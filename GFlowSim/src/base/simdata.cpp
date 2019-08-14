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
      if (v) dealloc_array_2d(v);
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

    // Get data from topology
    #if USE_MPI == 1
    neighbor_ranks = topology->get_neighbor_ranks();
    int n_size = neighbor_ranks.size();
    if (n_size>0) {
      // Send ids array.
      send_ids.resize(n_size);
      // Send ghost list array.
      send_ghost_list.resize(n_size);
      // Request list array.
      request_list.resize(n_size);
      // Recieve ghost sizes vector
      recv_ghost_sizes.resize(n_size);
      // Send size.
      send_size.resize(n_size);
      // Set up neighbor map
      for (int i=0; i<n_size; ++i) neighbor_map.insert(pair<int, int>(neighbor_ranks[i], i));
      // Set up buffer list.
      buffer_list.resize(n_size);
      // Set up recv_buffer.
      recv_buffer.resize(n_size);
    }
    // Set data width - position (sim_dimensions), velocity (sim_dimensions), radius (1), inverse mass (1), and type (1).
    data_width = 2*sim_dimensions + 3;
    // Set ghost data width - position (sim_dimensions).
    ghost_data_width = sim_dimensions;
    #endif

    //*****
    #if USE_MPI == 1
      if (topology->getNumProc()>1) {
        // Temporary fix - get rid of all particles not hosted on this processor.
        if (_size>0 && topology->is_initialized() ) {
          int rank = topology->getRank();
          int numProc = topology->getNumProc();
          for (int i=0; i<_size; ++i) {
            if (Type(i)>-1 && !topology->owned_particle(X(i))) markForRemoval(i);
          }
        }
      }
    #endif
    //*****

    // Sort the particles by position
    sortParticles();
  }

  void SimData::pre_integrate() {
    send_timer.clear(); 
    recv_timer.clear(); 
    barrier_timer.clear(); 
    ghost_send_timer.clear(); 
    ghost_recv_timer.clear(); 
    ghost_wait_timer.clear();
    exchange_search_timer.clear(); 
    ghost_search_timer.clear();
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
      if (v) dealloc_array_2d(v);
      v = alloc_array_2d<RealType>(num, sim_dimensions);
    }
    for (auto &s : sdata) {
      if (s) delete [] s;
      s = new RealType[num]; // Valgrind says there is an error here.
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

  void SimData::addParticle() {
    if (_size+1 > _capacity) {
      int new_capacity = max(32, static_cast<int>(0.25*_capacity));
      resize_owned(new_capacity);
    }
    // Reset all data
    reset_particle(_size);
    // Set type, give a global id
    Type(_size) = 0;
    id_map.insert(IPair(_size, next_global_id));
    Id(_size) = next_global_id++;
    ++_number;
    ++_size;

    // Assumes that there are no halo or ghost particles.
    ++_first_halo;
    ++_first_ghost;
  }

  //! \param num The number of particle slots to add.
  void SimData::addParticle(int num) {
    if (num<=0) return;
    if (_size+num > _capacity) {
      int new_capacity = max(32, static_cast<int>(0.25*(num+_size-_capacity)));
      resize_owned(new_capacity);
    }
    for (int i=0; i<num; ++i) {
      // Reset all data
      reset_particle(_size);
      // Set type, give a global id
      Type(_size) = 0;
      id_map.insert(IPair(_size, next_global_id));
      Id(_size) = next_global_id++;
      ++_number;
      ++_size;

      // Assumes that there are no halo or ghost particles.
      ++_first_halo;
      ++_first_ghost;
    }
  }

  //! \param x The position of the particle.
  //! \param v The velocity of the particle.
  //! \param sg The cutoff radius of the particle.
  //! \param im The inverse mass of the particle.
  //! \param type The type of the particle.
  void SimData::addParticle(const RealType *x, const RealType *v, const RealType sg, const RealType im, const int type) {
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
    id_map.insert(IPair(_size, next_global_id));
    Id(_size) = next_global_id++;
    ++_number;
    ++_size;

    // Assumes that there are no halo or ghost particles.
    ++_first_halo;
    ++_first_ghost;
  }

  void SimData::markForRemoval(const int id) {
    // If the particle has already been marked for removal, or is beyond the end of the array, return.
    if (Type(id)<0 || id>=_size) return;
    // Mark for removal, clear some data
    remove_list.insert(id);
    Type(id) = -1;
    id_map.erase(Id(id));
    Id(id) = -1;
    // This is probably unneccesary, but set V, F to zero.
    zeroVec(V(id), sim_dimensions);
    zeroVec(F(id), sim_dimensions);
  }

  void SimData::doParticleRemoval() {
    // If there is nothing to remove, we're done
    if (remove_list.empty() || _number==0) return;

    // Start simdata timer.
    timer.start();

    // Fill in all holes
    int count_back = _size, removed = 0;
    for(auto id : remove_list) {
      // The type of a particle that has been marked for removal is -1.
      do {
        --count_back;
      } while (Type(count_back)<0);

      if (count_back>id) {
        swap_particle(count_back, id); //** Formerly move_particle
        ++removed;
      }
      else break;
    }
    // Decrease number.
    _number -= remove_list.size();
    // Array is compressed.
    _size = _number;
    // Clear the remove list
    remove_list.clear();

    // We need to update. Removed could be zero if, e.g. only and all of the N particles were removed.
    // Arguably, you still might want to remake. There could be extra entries in the verlet list that 
    // it would be better to get rid of.
    if (removed>0) setNeedsRemake();

    // Start simdata timer.
    timer.stop();
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

  void SimData::updateHaloParticles() {
    // Start simdata timer.
    timer.start();

    // First pass: update the primary (actual) particle from the force data of the halo particles.
    // Doing this in two passes takes care of the fact that some particles may generate multiple halos.
    // We only have to update the forces, and then let the integrator take care of the rest.
    for (int i=0; i<halo_map.size(); i+=2) {
      int hid = halo_map[i]; // Halo id.
      int pid = halo_map[i+1]; // Primary id.
      // Update force
      plusEqVec(vdata[2][pid], vdata[2][hid], sim_dimensions);
    }
    // Second pass: update the force of the halo particles to match that of the primary particle.
    for (int i=0; i<halo_map.size(); i+=2) {
      int hid = halo_map[i]; // Halo id.
      int pid = halo_map[i+1]; // Primary id.
      // Update force
      copyVec(vdata[2][pid], vdata[2][hid], sim_dimensions);
    }

    // Start simdata timer.
    timer.stop();
  }

  void SimData::updateGhostParticles() {
  #if USE_MPI == 1
    // Return if the arrays are not allocated, or this there is only one processor.
    if (send_ids.empty() || topology->getNumProc()==1) return;
    // Update ghost particles.
    update_ghost_particles();
  #endif // USE_MPI == 1
  }

  RealType** SimData::X() {
    return vdata[0];
  }

  RealType* SimData::X_arr() {
    return *vdata[0];
  }

  RealType* SimData::X(int i) {
    return vdata[0][i];
  }

  RealType& SimData::X(int i, int d) {
    return vdata[0][i][d];
  }

  RealType** SimData::V() {
    return vdata[1];
  }

  RealType* SimData::V_arr() {
    return *vdata[1];
  }

  RealType* SimData::V(int i) {
    return vdata[1][i];
  }

  RealType& SimData::V(int i, int d) {
    return vdata[1][i][d];
  }

  RealType** SimData::F() {
    return vdata[2];
  }

  RealType* SimData::F_arr() {
    return *vdata[2];
  }

  RealType* SimData::F(int i) {
    return vdata[2][i];
  }

  RealType& SimData::F(int i, int d) {
    return vdata[2][i][d];
  }

  RealType** SimData::VectorData(int i) {
    if (i<vdata.size()) return vdata[i];
    else return nullptr;
  }

  RealType* SimData::Sg() {
    return sdata[0];
  }

  RealType& SimData::Sg(int i) {
    return sdata[0][i];
  }

  RealType* SimData::Im() {
    return sdata[1];
  }

  RealType& SimData::Im(int i) {
    return sdata[1][i];
  }

  RealType* SimData::ScalarData(int i) {
    if (i<sdata.size()) return sdata[i];
    else return nullptr;
  }

  int* SimData::Type() {
    return idata[0];
  }

  int& SimData::Type(int i) {
    return idata[0][i];
  }

  int* SimData::Id() {
    return idata[1];
  }

  int& SimData::Id(int i) {
    return idata[1][i];
  }

  int* SimData::IntegerData(int i) {
    if (i<idata.size()) return idata[i];
    else return nullptr;
  }

  const RealType** SimData::X() const {
    return const_cast<const RealType**>(vdata[0]);
  }

  const RealType* SimData::X_arr() const {
    return *vdata[0];
  }

  const RealType* SimData::X(int i) const {
    return vdata[0][i];
  }

  const RealType& SimData::X(int i, int d) const {
    return vdata[0][i][d];
  }

  const RealType** SimData::V() const {
    return const_cast<const RealType**>(vdata[1]);
  }

  const RealType* SimData::V_arr() const {
    return *vdata[1];
  }

  const RealType* SimData::V(int i) const {
    return vdata[1][i];
  }

  const RealType& SimData::V(int i, int d) const {
    return vdata[1][i][d];
  }

  const RealType** SimData::F() const {
    return const_cast<const RealType**>(vdata[2]);
  }

  const RealType* SimData::F_arr() const {
    return *vdata[2];
  }

  const RealType* SimData::F(int i) const {
    return vdata[2][i];
  }

  const RealType& SimData::F(int i, int d) const {
    return vdata[2][i][d];
  }

  const RealType* SimData::Sg() const {
    return sdata[0];
  }

  const RealType& SimData::Sg(int i) const {
    return sdata[0][i];
  }

  const RealType* SimData::Im() const {
    return sdata[1];
  }

  const RealType& SimData::Im(int i) const {
    return sdata[1][i];
  }

  const int* SimData::Type() const {
    return idata[0];
  }

  const int& SimData::Type(int i) const {
    return idata[0][i];
  }

  const int* SimData::Id() const {
    return idata[1];
  }

  const int& SimData::Id(int i) const {
    return idata[1][i];
  }

  bool SimData::Valid(int i) const {
    return -1<idata[1][i];
  }

  int SimData::requestVectorData(string name) {
    // Check if the data already exists
    auto it = vector_data_map.find(name);
    if (it!=vector_data_map.end()) return it->second;
    // Otherwise, create a data entry
    vector_data_map.insert(SIPair(name, vdata.size()));
    RealType **address = alloc_array_2d<RealType>(_capacity, sim_dimensions);
    vdata.push_back(address);
    // Return the entry
    return vdata.size()-1;
  }

  int SimData::requestScalarData(string name) {
    // Check if the data already exists
    auto it = scalar_data_map.find(name);
    if (it!=scalar_data_map.end()) return it->second;
    // Otherwise, create a data entry
    scalar_data_map.insert(SIPair(name, sdata.size()));
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
    integer_data_map.insert(SIPair(name, idata.size()));
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
    Vec x(sim_dimensions); x = X(id);
    Vec v(sim_dimensions); v = V(id);
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
    timer.start();

    // Mark halo particles for removal
    for (int i=0; i<halo_map.size(); i+=2) 
      markForRemoval(halo_map[i]);
    // Clear halo data
    halo_map.clear();
    // Clear the halo particle counter.
    _number_halo = 0;

    // Stop timer.
    timer.stop();
  }

  void SimData::removeGhostParticles() {
  #if USE_MPI == 1
    // Start timer.
    timer.start();

    // Check if this is necessary.
    if (n_ghosts<=0) return;
    // Remove ghost particles.
    for (int i=0; i<n_ghosts; ++i) markForRemoval(_first_ghost+i);

    // Stop timer.
    timer.start();
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
    Base::gflow->wrapPositions();

    #if USE_MPI == 1
      // Get rank and number of processors.
      int rank = topology->getRank();
      int numProc = topology->getNumProc();

      #if _CLANG_ == 1
        // If there are multiple processors.
        if (numProc>1 && !neighbor_ranks.empty()) {
          
          // --- Move particles that belong to other domains to those domains, and delete them from here. Then receive
          //     particles from other domains that belong to this domain.
          exchange_particles();

          // Start mpi timer.
          gflow->startMPIExchangeTimer();
          // Make sure all requests are processed.
          MPIObject::barrier(barrier_timer);
          // Stop mpi timer.
          gflow->stopMPIExchangeTimer();

          // Start adding ghost particles here.
          int save_first_ghost = _size;

          // --- Look for particles that need to be ghosts on other processors.
          if (gflow->use_ghosts()) create_ghost_particles();

          // Adding new particles increments _first_halo and _first_ghost, so we have to correct them.
          _first_ghost = save_first_ghost;
          _first_halo  = save_first_ghost;
        }
        // If there is only one processor, even though this is an MPI run.
        else doParticleRemoval();
      #else
        // \todo Cover this case later (non CLANG). Best way to do this is make objects that work with or without CLANG.
      #endif // _CLANG_ == 1
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

  int SimData::ntypes() const {
    return _ntypes;
  }

  void SimData::clearV() {
    RealType *v = V_arr();
    for (int i=0; i<_size*sim_dimensions; ++i) v[i] = 0;
  }

  void SimData::clearF() {
    RealType *f = F_arr();
    for (int i=0; i<_size*sim_dimensions; ++i) f[i] = 0;
  }

  void SimData::clearScalar(const string id) {
    int address = getScalarData(id);
    if (address>=0) {
      RealType *entry = ScalarData(address);
      for (int i=0; i<_size; ++i) entry[i] = 0;
    }
  }

  int SimData::getLocalID(int global) const {
    auto it = id_map.find(global);
    // Return the global iterator. We use -1 to mean "no such particle."
    return it==id_map.end() ? -1 : it->second;
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

  int SimData::getLastNGhostsSent() {
    return _last_n_ghosts_sent;
  }

  int SimData::getLastNGhostsRecv() {
    return _last_n_ghosts_recv;
  }

  int SimData::getLastNExchangeSent() {
    return _last_n_exchange_sent;
  }

  int SimData::getLastNExchangeRecv() {
    return _last_n_exchange_recv;
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
    vector_data_map.insert(SIPair(name, vector_data_map.size()));
    vdata.push_back(nullptr);
  }

  void SimData::addScalarData(string name) {
    scalar_data_map.insert(SIPair(name, scalar_data_map.size()));
    sdata.push_back(nullptr);
  }

  void SimData::addIntegerData(string name) {
    integer_data_map.insert(SIPair(name, integer_data_map.size()));
    idata.push_back(nullptr);
  }

  void SimData::resize_owned(int num) {
    // Compute new capacity
    int new_capacity = _capacity + num;
    // Allocate new vector data arrays
    for (auto &v : vdata) {
      RealType **nv = alloc_array_2d<RealType>(new_capacity, sim_dimensions);
      // Delete old array, set new
      if (v) {
  	    // Transfer data
  	    copyVec(*v, *nv, _size*sim_dimensions);
        // Delete old
        dealloc_array_2d(v);
      }
      // Initialize the reset of the data
      setVec(*nv, _size*sim_dimensions, new_capacity*sim_dimensions, static_cast<RealType>(0.));
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
    for (auto v : vdata) zeroVec(v[id], sim_dimensions);
    for (auto s : sdata) s[id] = 0.;
    for (auto i : idata) i[id] = -1;
  }

  void SimData::swap_particle(int id1, int id2) {
    // Get global IDs
    int g1 = Id(id1);
    int g2 = Id(id2);

    // Transfer data
    for (auto v : vdata) swapVec(v[id1], v[id2], sim_dimensions);
    for (auto s : sdata) std::swap(s[id1], s[id2]);
    for (auto i : idata) std::swap(i[id1], i[id2]);
    
    // Swap global ids
    auto it1 = id_map.find(g1);
    auto it2 = id_map.find(g2);
    if (id_map.end()!=it1 && id_map.end()!=it2) 
      std::swap(it1->second, it2->second);

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

#if USE_MPI == 1

  inline void SimData::exchange_particles() {
    // Start mpi timer.
    gflow->startMPIExchangeTimer();

    // Reset counters
    _last_n_exchange_sent = _last_n_exchange_recv = 0;

    // Get the bounds this processor manages.
    Bounds bounds = topology->getProcessBounds();

    // Clear send_ids buffer
    for (int i=0; i<neighbor_ranks.size(); ++i) send_ids[i].clear();
    
    /// --- Go through all particles, deciding which ones should migrate to other processors.
    exchange_search_timer.start();
    for (int id=0; id < _size; ++id) {
      if (Valid(id) && !bounds.contains(X(id))) {
        // Check which processor the particle actually belongs on. We can use the send_ids buffer since we are going to 
        // clear it anyways.
        int n_rank = topology->domain_ownership(X(id));
        // Find which entry the id should be put into to go to the correct neighbor.
        auto it = neighbor_map.find(n_rank);
        // If there are repulsive boundaries, extend the simulation bounds for the node in those directions, so particles don't go "out of bounds"
        // when they pass over the boundary a little to get repulsed back in. In this case, X can evaluate to not being within the bounds, even
        // though it should stay within these bounds.
        if (it!=neighbor_map.end()) {
          send_ids[it->second].push_back(id);
          ++_last_n_exchange_sent;
        }
      }
    }
    exchange_search_timer.stop();

    // Send particle information, deleting the particles that we send.
    send_timer.start();
    for (int i=0; i<neighbor_ranks.size(); ++i)
      send_particle_data(send_ids[i], neighbor_ranks[i], buffer_list[i], &request_list[i], &request, true);
    send_timer.stop();

    // Stop mpi timer. Particle removal counts as a simdata update task.
    gflow->stopMPIExchangeTimer();

    // Do particle removal. We do this here so we get rid of all the particles that have migrated to other processors. There will be no "holes"
    // in the array after the particle removal.
    doParticleRemoval(); 

    // Stop mpi timer.
    gflow->startMPIExchangeTimer();
    
    // --- Receive particles from other processors.

    // Recieve particle information, and use it to create new particles.
    recv_timer.start();
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      int n_rank = neighbor_ranks[i];
      _last_n_exchange_recv += recv_new_particle_data(n_rank, recv_buffer[i]);
    }
    recv_timer.stop();

    // Stop mpi timer.
    gflow->stopMPIExchangeTimer();
  }

  inline void SimData::create_ghost_particles() {
    // Start mpi timer.
    gflow->startMPIGhostTimer();

    // First, clear all the send_ghost_lists entries.
    for (int i=0; i<neighbor_ranks.size(); ++i) send_ghost_list[i].clear();

    ghost_search_timer.start();
    vector<int> overlaps; // Helping vector
    RealType skin_depth = handler->getSkinDepth();
    for (int id=0; id < _size; ++id) {
      // The particle cutoff that should be used to test whether the particle is close enough to another domain.
      // \todo There may be a better/more correct way to do this.
      RealType cutoff = 2 * Sg(id) * forceMaster->getMaxCutoff(Type(id)) + skin_depth;
      // Check if the particle overlaps with another domain.
      topology->domain_overlaps(X(id), cutoff, overlaps);

      // Store the particle id in the send_ghost_list entry for every processor we need to send this particle to as a ghost.
      for (auto proc_n : overlaps) send_ghost_list[proc_n].push_back(id);            
    }
    ghost_search_timer.stop();

    // --- Make sure the entries of buffer list are large enough.
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      // Find the number of elements we need to send. The total amount of space we need will be s*ghost_data_width, since each
      // particle needs size ghost_data_width.
      int s = send_ids[i].size();
      // If the i-th buffer is smaller than it needs to be, resize it.
      if (buffer_list[i].size()<s*ghost_data_width) buffer_list[i].resize(s*ghost_data_width);
    }

    // --- Tell neighboring processors how many particles to expect.
    ghost_send_timer.start();
    // Send particle information, but do not delete the original particles, since they will be ghosts on the other processors.
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      send_particle_data(send_ghost_list[i], neighbor_ranks[i], buffer_list[i], &request_list[i], &request, false);
    }
    ghost_send_timer.stop();

    // Reset n_ghosts.
    n_ghosts = 0;
    // Get ghosts from all neighbors.
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      // Rank of the i-th neighbor.
      int n_rank = neighbor_ranks[i];

      // Recieve ghost particles, create new particles for them.
      ghost_recv_timer.start();
      recv_ghost_sizes[i] = recv_new_particle_data(n_rank, recv_buffer[i]);
      ghost_recv_timer.stop();
      
      // Update number of ghosts.
      n_ghosts += recv_ghost_sizes[i];
    }

    // Stop mpi timer.
    gflow->stopMPIGhostTimer();
  }

  inline void SimData::update_ghost_particles() {
    // Start mpi timer.
    gflow->startMPIGhostTimer();

    // Reset counters.
    _last_n_ghosts_sent = _last_n_ghosts_recv = 0;

    // Update the positions information of ghost particles on other processors.
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      // How many ghost particles are hosted on the i-th processor, and should have their information returned to there.
      int size = send_ghost_list[i].size();
      // Only send data if there is data to send.
      if (size>0) {
        // Update counter.
        _last_n_ghosts_sent += size;
        // Make sure buffer is large enough.
        if (buffer_list[i].size()<size*ghost_data_width) buffer_list[i].resize(size*ghost_data_width);
        // Pack the data into buffer_list
        int j = 0;
        for (int id : send_ghost_list[i]) {
          copyVec(X(id), &buffer_list[i][ghost_data_width*j], sim_dimensions);
          ++j;
        }
        
        // Send the data (non-blocking) back to the processor.
        ghost_send_timer.start();
        MPI_Isend(buffer_list[i].data(), size*ghost_data_width, MPI_FLOAT, neighbor_ranks[i], 0, MPI_COMM_WORLD, &request);
        ghost_send_timer.stop();
      }
    }

    // Start non-blocking receives of ghost particle data.
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      int size = recv_ghost_sizes[i];
      // Only expect a message if size is positive.
      if (size>0) {
        _last_n_ghosts_recv += size;
        // Make sure buffer size is good.
        if (recv_buffer[i].size() < size*ghost_data_width) recv_buffer[i].resize(size*ghost_data_width);
        // Start the non-blocking receieve.
        ghost_recv_timer.start();
        MPI_Irecv(recv_buffer[i].data(), size*ghost_data_width, MPI_FLOAT, neighbor_ranks[i], 0, MPI_COMM_WORLD, &request_list[i]);
        ghost_recv_timer.stop();
      }
    }

    // Collect received data and pack it.
    for (int i=0, p_id=_first_ghost; i<neighbor_ranks.size(); ++i) {
      int size = recv_ghost_sizes[i];
      if (size>0) {



        // Wait for request to be filled.
        MPIObject::wait(request_list[i], ghost_wait_timer);

        // ghost_wait_timer.start();
        // int count = 0;
        // while (!MPIObject::test(request_list[i])) ++count;
        // ghost_wait_timer.stop();
        // cout << "Rank: " << topology->getRank() << ", Iter: " << gflow->getIter() << ", Flag: " << gflow->handler_remade() << ", " << count << "\n";

        // Unpack the position data into the ghost particles.
        for (int j=0; j<size; ++j, ++p_id) {
          copyVec(&recv_buffer[i][ghost_data_width*j], X(p_id), sim_dimensions);
        }
      }
    }

    // Stop mpi timer.
    gflow->stopMPIGhostTimer();
  }

  inline void SimData::send_particle_data(const vector<int>& send_ids, int n_rank, vector<RealType>& buffer, MPI_Request* size_request, MPI_Request* send_request, bool remove) {
    // Tell the processor how many particles to expect. Use a non-blocking send, store the request 
    int size = send_ids.size();
    // MPIObject::send_single(size, n_rank);
    MPI_Isend(&size, 1, MPI_INT, n_rank, 0, MPI_COMM_WORLD, size_request);    
    // Send the actual particles, if there are any.
    if (size>0) {
      // Make sure buffer is big enough to send data.
      if (buffer.size()<size*data_width) buffer.resize(size*data_width);
      // Send the actual data. Copy data into buffer
      for (int j=0; j<size; ++j) {
        int id = send_ids[j];
        // Copy particle information to the buffer. \todo Automate a way to specify arbitrary subsets of the particle data to send.
        copyVec(X(id), &buffer[data_width*j], sim_dimensions); // Position
        copyVec(V(id), &buffer[data_width*j + sim_dimensions], sim_dimensions); // Velocity
        buffer[data_width*j + 2*sim_dimensions + 0] = Sg(id); // Radius
        buffer[data_width*j + 2*sim_dimensions + 1] = Im(id); // Inverse mass
        buffer[data_width*j + 2*sim_dimensions + 2] = Type(id); // Type
        // Mark particle for removal.
        if (remove) markForRemoval(id);
      }
      // Send the data (non-blocking).
      MPI_Isend(buffer.data(), size*data_width, MPI_FLOAT, n_rank, 0, MPI_COMM_WORLD, send_request);
    }
  }

  inline int SimData::recv_new_particle_data(int n_rank, vector<RealType>& buffer) {
    // Tell the processor how many particles to expect
    int size = 0;
    // MPIObject::recv_single(size, n_rank);
    MPI_Recv(&size, 1, MPI_INT, n_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);     
    // Get the actual particles, if there are any.
    if (size>0) {
      // Vectors for unpacking.
      Vec X(sim_dimensions), V(sim_dimensions);
      // Resize buffer if neccessary.
      if (buffer.size()<size*data_width) buffer.resize(size*data_width); 
      // Receive buffer.
      MPI_Recv(buffer.data(), size*data_width, MPI_FLOAT, n_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Add particle.
      for (int j=0; j<size; ++j) {
        // Position, (velocity is zero), 
        copyVec(&buffer[data_width*j], X);
        copyVec(&buffer[data_width*j + sim_dimensions], V);
        RealType r    = buffer[data_width*j + 2*sim_dimensions + 0];
        RealType im   = buffer[data_width*j + 2*sim_dimensions + 1];
        RealType type = buffer[data_width*j + 2*sim_dimensions + 2];

        // CHECK
        int t = static_cast<int>(type);
        if (t>0) throw false;

        // Add particle to simdata.
        addParticle(X.data, V.data, r, im, static_cast<int>(type));
      } 
    }
    // Return the number of particles that were received.
    return size;
  }

#endif // USE_MPI == 1

}
