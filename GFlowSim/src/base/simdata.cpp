#include "simdata.hpp"
// Other files
#include "integrator.hpp"
#include "forcemaster.hpp"
#include "datamaster.hpp"
#include "interactionhandler.hpp"
#include "topology.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  SimData::SimData(GFlow *gflow) : Base(gflow) {
    // Add default vector data entries
    addVectorData("X");
    addVectorData("V");
    addVectorData("F");
    // Add default scalar data entries
    addScalarData("Sg");
    addScalarData("Im");
    // Add default integer data entries
    addIntegerData("Type");
    addIntegerData("ID");
  }

  SimData::~SimData() {
    for (auto &v : vdata) 
      if (v) delete [] v;
    vdata.clear();
    for (auto &s : sdata)
      if (s) delete [] s;
    sdata.clear();
    for (auto &i : idata)
      if (i) delete [] i;
    idata.clear();
  }

  //! @brief Initialize the atom container.
  void SimData::initialize() {
    // Call base's initialize
    Base::initialize();

    // Get ntypes from force master
    _ntypes = forceMaster->getNTypes();

    // Set data width - send all data to adjacent processors when migrating ghosts.
    // \todo We don't actually need to migrate force data.
    data_width = vdata.size()*sim_dimensions + sdata.size() + idata.size();

    // Set ghost data width - position (sim_dimensions).
    ghost_data_width = sim_dimensions;
    // Do we need to send velocity?
    if (send_ghost_velocity) ghost_data_width += sim_dimensions;
    // Do we need to send angular velocity?
    if (send_ghost_omega) ghost_data_width += 1; // In two dimensions, omega is a scalar.

    #if USE_MPI == 1
    if (topology->getNumProc()>1) {
      // Temporary fix - get rid of all particles not hosted on this processor.
      if (_size>0 && topology->is_initialized() ) {
        int rank = topology->getRank();
        int numProc = topology->getNumProc();
        for (int i=0; i<_size; ++i) {
          if (Type(i)>-1 && !topology->owned_particle(X(i)))
            markForRemoval(i);
        }
      }
    }
    #endif

    // Sort the particles by position
    sortParticles();
  }

  void SimData::pre_integrate() {
    // Clear the timed object timer.
    clear_timer();
    // Clear other timers.
    send_timer.clear_timer(); 
    recv_timer.clear_timer(); 
    barrier_timer.clear_timer(); 
    ghost_send_timer.clear_timer(); 
    ghost_recv_timer.clear_timer(); 
    ghost_wait_timer.clear_timer();
    exchange_search_timer.clear_timer(); 
    ghost_search_timer.clear_timer();
    // Remove any bad particles.
    removeBadParticles();

    // Make sure data widths are up to date.
    #if USE_MPI == 1
    // Set data width - send all data to adjacent processors when migrating ghosts.
    // \todo We don't actually need to migrate force data.
    data_width = vdata.size()*sim_dimensions + sdata.size() + idata.size();
    // Set ghost data width - position (sim_dimensions).
    ghost_data_width = sim_dimensions;
    // Do we need to send velocity?
    if (send_ghost_velocity) ghost_data_width += sim_dimensions;
    // Do we need to send angular velocity?
    if (send_ghost_omega) ghost_data_width += 1; // In two dimensions, omega is a scalar.
    #endif // USE_MPI == 1
  }

  void SimData::post_integrate() {
    // Mark extraneous particles for removal
    removeHaloAndGhostParticles();

    // Do actual particle removal.
    doParticleRemoval();

    // Set markers.
    _first_halo = _first_ghost = _size;    
  }

  //! @brief Reserve space for particles, extending the lengths of all arrays to the requested size.
  void SimData::reserve(int num) {
    for (auto &v : vdata) {
      if (v) delete [] v;
      v = new real[num*sim_dimensions];
    }
    for (auto &s : sdata) {
      if (s) delete [] s;
      s = new real[num]; // Valgrind says there is an error here.
    }
    for (auto &i : idata) {
      if (i) delete [] i;
      i = new int[num];
    }
    // Reset numbers
    _number = 0;
    _size = 0;
    _capacity = num;
  }

  int SimData::addParticle() {
    if (_size+1 > _capacity) {
      int new_capacity = max(32, static_cast<int>(0.25*_capacity));
      resize_owned(new_capacity);
    }
    // Reset all data
    reset_particle(_size);
    // Set type, give a global id
    Type(_size) = 0;
    if (use_id_map) id_map.emplace(_size, next_global_id);
    Id(_size) = next_global_id++;
    ++_number;
    ++_size;

    // Assumes that there are no halo or ghost particles.
    ++_first_halo;
    ++_first_ghost;
    // Return particle id.
    return _size-1;
  }

  //! \param num The number of particle slots to add.
  int SimData::addParticle(int num) {
    if (num<=0) return -1;
    if (_size+num > _capacity) {
      int new_capacity = max(32, static_cast<int>(0.25*(num+_size-_capacity)));
      resize_owned(new_capacity);
    }
    int address = _size;
    for (int i=0; i<num; ++i) {
      // Reset all data
      reset_particle(_size);
      // Set type, give a global id
      Type(_size) = 0;
      if (use_id_map) id_map.emplace(_size, next_global_id);
      Id(_size) = next_global_id++;
      ++_number;
      ++_size;

      // Assumes that there are no halo or ghost particles.
      ++_first_halo;
      ++_first_ghost;
    }
    // Return first particle id.
    return address;
  }

  //! \param x The position of the particle.
  //! \param v The velocity of the particle.
  //! \param sg The cutoff radius of the particle.
  //! \param im The inverse mass of the particle.
  //! \param type The type of the particle.
  int SimData::addParticle(const RealType *x, const RealType *v, const RealType sg, const RealType im, const int type) {
    // If not enough spots to add a new owned particle, create more
    if (_size+1 > _capacity) {
      int new_capacity = max(32, static_cast<int>(0.25*_size));
      resize_owned(new_capacity);
    }
    // Reset all data
    reset_particle(_size);
    // Set data
    copyVec(x, X(_size), sim_dimensions);
    copyVec(v, V(_size), sim_dimensions);
    Sg(_size) = sg;
    Im(_size) = im;
    Type(_size) = type;
    if (use_id_map) id_map.emplace(_size, next_global_id);
    Id(_size) = next_global_id++;
    ++_number;
    ++_size;

    // Assumes that there are no halo or ghost particles.
    ++_first_halo;
    ++_first_ghost;
    // Return particle id.
    return _size-1;
  }

  void SimData::markForRemoval(const int id) {
    // If the particle has already been marked for removal, or is beyond the end of the array, return.
    if (Type(id)<0 || id>=_size) return;
    // Mark for removal, clear some data
    remove_list.emplace(id);
    Type(id) = -1;
    if (use_id_map) id_map.erase(Id(id));
    Id(id) = -1;
    // This is probably unneccesary, but set V, F to zero.
    zeroVec(V(id), sim_dimensions);
    zeroVec(F(id), sim_dimensions);
  }

  void SimData::doParticleRemoval() {
    // If there is nothing to remove, we're done
    if (remove_list.empty() || _number==0) return;

    // Start simdata timer.
    start_timer();

    // Fill in all holes
    int count_back = _size, removed = 0;
    for(auto id : remove_list) {
      // The type of a particle that has been marked for removal is -1.
      do {
        --count_back;
      } while (Type(count_back)<0);

      if (count_back>id) {
        swap_particle(count_back, id);
        ++removed;
      }
      else break;
    }
    // Decrease number.
    _number -= remove_list.size();
    // Array is compressed.
    _size = _number;
    _first_ghost = _size;
    _first_halo = _size;
    // Clear the remove list
    remove_list.clear();

    // We need to update. Removed could be zero if, e.g. only and all of the N particles were removed.
    // Arguably, you still might want to remake. There could be extra entries in the verlet list that 
    // it would be better to get rid of.
    if (removed>0) setNeedsRemake();

    // Start simdata timer.
    stop_timer();
  }

  void SimData::sortParticles() {
    // We must first remove halo and ghost particles.
    removeHaloAndGhostParticles();

    // Make sure all particles are valid, and compressed
    doParticleRemoval(); // This only sets the needs remake flag if it removes particles.

    // Quick sort
    quick_sort(0, _number-1, 0);
    recursion_help (0, _number-1, 1);

    // Set needs remake flag
    setNeedsRemake(true);
  }

  void SimData::sortParticles(Vec& direction) {
    // Make sure all particles are valid, and compressed
    doParticleRemoval(); // This only sets the needs remake flag if it removes particles.
    // FOR NOW: JUST SORT ALONG X AXIS
    quick_sort(0, _number-1, 0);
    // Set needs remake flag
    setNeedsRemake(true);
  }

  bool SimData::removeBadParticles(bool do_removal) {
    // Most of the time, we probably won't need to remove any particles. We only need to remake structures if we do. So keep track.
    bool removed_some = false;
    // Look for bad particles.
    for (int i=0; i<_size; ++i) {
      for (int d=0; d<sim_dimensions; ++d) {
        // If a component of position or velocity is nan, remove that particle.
        if (isnan(X(i, d)) || isnan(V(i, d))) {
          markForRemoval(i);
          removed_some = true;
          break;
        }
      }
    }
    if (removed_some && do_removal) {
      // Remove particles.
      doParticleRemoval();
      // Set the needs remake flag
      setNeedsRemake(true);
    }
    return removed_some;
  }

  void SimData::updateHaloParticles() {
    // Start simdata timer.
    start_timer();

    // First pass: update the primary (actual) particle from the force data of the halo particles.
    // Doing this in two passes takes care of the fact that some particles may generate multiple halos.
    // We only have to update the forces, and then let the integrator take care of the rest.
    auto f = F();
    for (int i=0; i<halo_map.size(); i+=2) {
      int hid = halo_map[i]; // Halo id.
      int pid = halo_map[i+1]; // Primary id.
      // Update force
      plusEqVec(f(pid), f(hid), sim_dimensions);
    }
    // Second pass: update the force of the halo particles to match that of the primary particle.
    for (int i=0; i<halo_map.size(); i+=2) {
      int hid = halo_map[i]; // Halo id.
      int pid = halo_map[i+1]; // Primary id.
      // Update force
      copyVec(f(pid), f(hid), sim_dimensions);
    }

    // Start simdata timer.
    stop_timer();
  }

  vec_access SimData::X() {
    return vec_access(vdata[0], sim_dimensions);
  }

  RealType* SimData::X(const int i) {
    return X()(i);
  }

  RealType& SimData::X(const int i, const int d) {
    return X()(i, d);
  }

  vec_access SimData::V() {
    return vec_access(vdata[1], sim_dimensions);
  }

  RealType* SimData::V(const int i) {
    return V()(i);
  }

  RealType& SimData::V(const int i, const int d) {
    return V()(i, d);
  }

  vec_access SimData::F() {
    return vec_access(vdata[2], sim_dimensions);
  }

  RealType* SimData::F(const int i) {
    return F()(i);
  }

  RealType& SimData::F(const int i, const int d) {
    return F()(i, d);
  }

  vec_access SimData::VectorData(const int i) {
    if (i<0 || vdata.size()<=i) throw false; // \todo Real error
    return vec_access(vdata[i], sim_dimensions);
  }

  vec_access SimData::VectorData(const string& name) {
    int i = getVectorData(name);
    return VectorData(i);
  }

  real* SimData::VectorData(const int i, const int id) {
    return VectorData(i)(id);
  }

  scalar_access SimData::Sg() {
    return scalar_access(sdata[0]);
  }

  RealType& SimData::Sg(const int i) {
    return Sg()(i);
  }

  scalar_access SimData::Im() {
    return scalar_access(sdata[1]);
  }

  RealType& SimData::Im(const int i) {
    return Im()(i);
  }

  scalar_access SimData::ScalarData(const int i) {
    // If i is bad, return a null access.
    if (i<0 || sdata.size()<=i) return scalar_access();
    return scalar_access(sdata[i]);
  }

  scalar_access SimData::ScalarData(const string& name) {
    int i = getScalarData(name);
    // Return access to the data.
    if (0<i) return ScalarData(i);
    // Return a null scalar access.
    else return scalar_access();
  }

  real& SimData::ScalarData(const int i, const int id) {
    return ScalarData(i)(id);
  }

  integer_access SimData::Type() {
    return integer_access(idata[0]);
  }

  int& SimData::Type(int i) {
    return Type()(i);
  }

  integer_access SimData::Id() {
    return integer_access(idata[1]);
  }

  int& SimData::Id(const int i) {
    return Id()(i);
  }

  integer_access SimData::IntegerData(const int i) {
    if (i<0 || idata.size()<=i) throw false; // \todo Real error
    return integer_access(idata[i]);
  }

  integer_access SimData::IntegerData(const string& name) {
    int i = getIntegerData(name);
    return IntegerData(i);
  }

  int& SimData::IntegerData(const int i, const int id) {
    return IntegerData(i)(id);
  }

  int SimData::requestVectorData(string name) {
    // Check if the data already exists
    auto it = vector_data_map.find(name);
    if (it!=vector_data_map.end()) return it->second;
    // Otherwise, create a data entry
    vector_data_map.emplace(name, vdata.size());
    vdata.push_back(new real[_capacity*sim_dimensions]);
    // Return the entry
    return vdata.size()-1;
  }

  int SimData::requestScalarData(string name) {
    // Check if the data already exists
    auto it = scalar_data_map.find(name);
    if (it!=scalar_data_map.end()) return it->second;
    // Otherwise, create a data entry
    scalar_data_map.emplace(name, sdata.size());
    RealType *address = new RealType[_capacity];
    sdata.push_back(address);
    // Return the entry
    return sdata.size()-1;
  }

  int SimData::requestIntegerData(string name) {
    // Check if the data already exists
    auto it = integer_data_map.find(name);
    if (it!=integer_data_map.end()) return it->second;
    // Otherwise, create a data entry
    integer_data_map.emplace(name, idata.size());
    int *address = new int[_capacity];
    idata.push_back(address);
    // Return the entry
    return idata.size()-1;
  }

  int SimData::getVectorData(string name) {
    // Check if the data already exists
    auto it = vector_data_map.find(name);
    if (it!=vector_data_map.end()) return it->second;
    else return -1;
  }

  int SimData::getScalarData(string name) {
    // Check if the data already exists
    auto it = scalar_data_map.find(name);
    if (it!=scalar_data_map.end()) return it->second;
    else return -1;
  }

  int SimData::getIntegerData(string name) {
    // Check if the data already exists
    auto it = integer_data_map.find(name);
    if (it!=integer_data_map.end()) return it->second;
    else return -1;
  }

  //! \param id The id of the particle to copy
  //! \param displacement How should the halo particle be displaced relative to the original particle.
  void SimData::createHaloOf(int id, const Vec& displacement) {
    // Record local ids
    halo_map.push_back(_size); // Halo local id
    halo_map.push_back(id);    // Original local id
    // Make a copy of the particle
    Vec x(sim_dimensions, X(id));
    Vec v(sim_dimensions, V(id));
    RealType radius = Sg(id);
    RealType im = Im(id);
    int type = Type(id);
    // Add a particle to be the halo particle.
    addParticle(x.data, v.data, radius, im, type);
    int id2 = _size-1;
    plusEqVec(X(id2), displacement.data, sim_dimensions); // Displace halo particle
    // Increment the halo particles counter.
    ++_number_halo;
    //! \todo There may be other data we should copy
  }

  void SimData::removeHaloParticles() {
    // Start timer.
    start_timer();

    // Mark halo particles for removal
    for (int i=0; i<halo_map.size(); i+=2) 
      markForRemoval(halo_map[i]);
    // Clear halo data
    halo_map.clear();
    // Clear the halo particle counter.
    _number_halo = 0;

    // Stop timer.
    stop_timer();
  }

  void SimData::removeGhostParticles() {
  #if USE_MPI == 1
    // Start timer.
    start_timer();

    // Check if this is necessary.
    if (_number_ghost<=0) return;
    // Remove ghost particles.
    for (int i=0; i<_number_ghost; ++i) markForRemoval(_first_ghost+i);

    // Stop timer.
    stop_timer();
  #endif
  }

  void SimData::removeHaloAndGhostParticles() {
    removeGhostParticles();
    removeHaloParticles();
  }

  void SimData::update() {
    /// Remove old halo and ghost particles
    removeHaloAndGhostParticles();

    // Wrap the particles, so they are in their cannonical positions
    gflow->wrapPositions();

    #if USE_MPI == 1
      // Get rank and number of processors.
      int rank = topology->getRank();
      int numProc = topology->getNumProc();

      // Sync the terminate flag.
      gflow->syncRunning();

      // If there are multiple processors.
      if (numProc>1) {

        // Move particles that belong to other domains to those domains, and delete them from here. Then receive
        // particles from other domains that belong to this domain.
        topology->exchange_particles();

        // Start adding ghost particles here.
        int save_first_ghost = _size;

        // --- Look for particles that need to be ghosts on other processors.
        if (gflow->use_ghosts()) {
          // Create ghost particles.
          topology->create_ghost_particles();
        }

        // Adding new particles increments _first_halo and _first_ghost, so we have to correct them.
        _first_ghost = save_first_ghost;
        _first_halo  = save_first_ghost;
      }
      // If there is only one processor, even though this is an MPI run.
      else doParticleRemoval();
    #else 
      // If not parallel, just do particle removal.
      doParticleRemoval();
    #endif // USE_MPI == 1
  }

  int SimData::size() const {
    return _size;
  }

  int SimData::size_owned() const {
    return _first_ghost;
  }

  int SimData::number() const {
    return _number;
  }

  int SimData::number_owned() const {
    return _number - _number_halo - _number_ghost;
  }

  int SimData::number_ghosts() const {
    return _number_ghost;
  }

  int SimData::first_ghost() const {
    return _first_ghost;
  }

  int SimData::ntypes() const {
    return _ntypes;
  }

  void SimData::clearV() {
    // Start simdata timer.
    start_timer();
    // Clear velocities
    auto v = V();
    for (int i=0; i<_size*sim_dimensions; ++i) v[i] = 0;
    // Stop simdata timer.
    stop_timer();
  }

  void SimData::clearF() {
    // Start simdata timer.
    start_timer();
    // Clear forces
    auto f = F();
    for (int i=0; i<_size*sim_dimensions; ++i) f[i] = 0;
    // Stop simdata timer.
    stop_timer();
  }

  void SimData::clearScalar(const string id) {
    // Start simdata timer.
    start_timer();
    // Clear scalar
    int address = getScalarData(id);
    if (address>=0) {
      auto entry = ScalarData(address);
      for (int i=0; i<_size; ++i) entry(i) = 0;
    }
    // Stop simdata timer.
    stop_timer();
  }

  int SimData::getLocalID(int global) const {
    if (use_id_map) {
      auto it = id_map.find(global);
      // Return the global iterator. We use -1 to mean "no such particle."
      return it==id_map.end() ? -1 : it->second;
    }
    return -2;
  }

  int SimData::getNextGlobalID() const {
    return next_global_id;
  }

  const BCFlag* SimData::getBCs() const {
    return Base::gflow->getBCs();
  }

  Bounds SimData::getBounds() const {
    return gflow->getBounds();
  }

  bool SimData::getNeedsRemake() {
    return needs_remake;
  }

  int SimData::getFirstHalo() {
    return _first_halo;
  }

  int SimData::getFirstGhost() {
    return _first_ghost;
  }

  bool SimData::isReal(int id) {
    return -1<id && id<_first_halo;
  }

  bool SimData::isHalo(int id) {
    return _first_halo <= id && id < _first_ghost;
  }

  bool SimData::isGhost(int id) {
    return _first_ghost <= id && id < _size;
  }

  void SimData::setFirstHalo(int id) {
    _first_halo = id;
  }

  void SimData::setFirstGhost(int id) {
    _first_ghost = id;
  }

  void SimData::setNeedsRemake(bool r) {
    // Set the local flag. \todo Only use flags in gflow?
    needs_remake = r;
    // Set the flag in simdata.
    gflow->simdata_needs_remake() = r;
    // If set to true, tell datamaster that local ids may have changed.
    if (r) dataMaster->setLocalsChanged(r);
  }

  void SimData::addVectorData(string name) {
    vector_data_map.emplace(name, vector_data_map.size());
    vdata.push_back(nullptr);
  }

  void SimData::addScalarData(string name) {
    scalar_data_map.emplace(name, scalar_data_map.size());
    sdata.push_back(nullptr);
  }

  void SimData::addIntegerData(string name) {
    integer_data_map.emplace(name, integer_data_map.size());
    idata.push_back(nullptr);
  }

  void SimData::resize_owned(int num) {
    // Compute new capacity
    int new_capacity = _capacity + num;
    // Allocate new vector data arrays
    for (auto &v : vdata) {
      RealType *nv = new real[new_capacity*sim_dimensions];
      // Delete old array, set new
      if (v) {
  	    // Transfer data
  	    copyVec(v, nv, _size*sim_dimensions);
        // Delete old
        delete [] v;
      }
      // Initialize the reset of the data
      setVec(nv, _size*sim_dimensions, new_capacity*sim_dimensions, static_cast<real>(0.));
      // Set pointer
      v = nv;
    }
    // Allocate new scalar data arrays
    for (auto &s : sdata) {
      RealType *ns = new RealType[new_capacity];
      // Delete old array, set new
      if (s) {
      	// Transfer data
      	copyVec(s, ns, _size);
      	// Delete old 
      	delete [] s;
      }
      // Initialize the rest of the data
      setVec(ns, _size, new_capacity, static_cast<RealType>(0.));
      // Set pointer
      s = ns;
    }
    // Allocate new integer data
    for (auto &i : idata) {
      int *ni = new int[new_capacity];
      // Delete old array, set new
      if (i) {
      	// Transfer data
      	copyVec(i, ni, _size);
      	delete [] i;
      }
      // Initialize the rest of the data
      setVec(ni, _size, new_capacity, -1); // Set to -1 so type will be -1
      // Set pointer
      i = ni;
    }
    // Set new sizes
    _capacity += num;
  }

  void SimData::reset_particle(int id) {
    for (int i=0; i<vdata.size(); ++i) zeroVec(VectorData(i)(id), sim_dimensions);
    for (auto s : sdata) s[id] = 0.;
    for (auto i : idata) i[id] = -1;
  }

  void SimData::swap_particle(int id1, int id2) {
    // Get global IDs
    int g1 = Id(id1);
    int g2 = Id(id2);

    // Transfer data
    for (int i=0; i<vdata.size(); ++i) {
      auto v = VectorData(i);
      swapVec(v(id1), v(id2), sim_dimensions);
    }
    for (int i=0; i<sdata.size(); ++i) {
      auto s = ScalarData(i);
      std::swap(s(id1), s(id2));
    }
    for (int i=0; i<idata.size(); ++i) {
      auto id = IntegerData(i);
      std::swap(id(id1), id(id2));
    }
    
    // Swap global ids
    if (use_id_map) {
      auto it1 = id_map.find(g1);
      auto it2 = id_map.find(g2);
      if (it1!=id_map.end()) it1->second = id2;
      if (it2!=id_map.end()) it2->second = id1;
    }

    // Set flag
    setNeedsRemake(true);
  }

  void SimData::quick_sort(int start, int end, int dim) {
    if (start<end && dim<sim_dimensions) {
      int partition = quick_sort_partition(start, end, dim);
      quick_sort(start, partition, dim);
      quick_sort(partition+1, end, dim);
    }
  }

  int SimData::quick_sort_partition(int start, int end, int dim) {
    RealType pivot = X((start + end)/2, dim);
    int i = start-1, j = end+1;
    while (true) {
      do {
        ++i;
      } while (X(i, dim)<pivot);
      do {
        --j;
      } while (X(j, dim)>pivot);      
      // Termination condition
      if (j<=i) return j;
      // If not, swap particles i and j
      swap_particle(i, j);
    }
  }

  void SimData::recursion_help(int start, int end, int dim) {
    if (dim>=sim_dimensions) return;
    // Do next level of quicksorts, for the next dimension
    const int sort_bins = 5, min_particles = 10;
    const int ds = (end-start)/sort_bins;

    if (ds>min_particles) {
      for (int i=0; i<sort_bins; ++i) {
        quick_sort(i*ds, (i+1)*ds, dim);
        recursion_help (i*ds, (i+1)*ds, dim+1);
      }
      // Potentially, there is some left over
      quick_sort(sort_bins*ds, _number-1, dim);
      recursion_help (sort_bins*ds, _number-1, dim+1);
    }
  }

}
