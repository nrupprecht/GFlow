#include "simdata.hpp"
// Other files
#include "../utility/memory.hpp"
#include "../base/forcemaster.hpp"
#include "../base/datamaster.hpp"

#include "integrator.hpp"

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

    // Sort the particles by position
    sortParticles();
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
  }

  void SimData::exchangeParticles() {
    // @todo Implement
  }

  void SimData::sortParticles() {
    // We must first remove halo and ghost particles.
    removeHaloAndGhostParticles();
    // Make sure all particles are valid, and compressed
    doParticleRemoval(); // This only sets the needs remake flag if it removes particles.

    // Quick sort
    quick_sort_help(0, _number-1, 0);
    recursion_help (0, _number-1, 1);
    // Set needs remake flag
    setNeedsRemake(true);
  }

  void SimData::sortParticles(Vec& direction) {
    // Make sure all particles are valid, and compressed
    doParticleRemoval(); // This only sets the needs remake flag if it removes particles.

    // FOR NOW: JUST SORT ALONG X AXIS
    quick_sort_help(0, _number-1, 0);

    // Set needs remake flag
    setNeedsRemake(true);
  }

  void SimData::updateHaloParticles() {
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
    // Mark halo particles for removal
    for (int i=0; i<halo_map.size(); i+=2) 
      markForRemoval(halo_map[i]);
    // Clear halo data
    halo_map.clear();
    // Clear the halo particle counter.
    _number_halo = 0;
  }

  void SimData::removeGhostParticles() {
    // Mark ghost particles for removal
    for (int i=0; i<ghost_map.size(); i+=2) 
      markForRemoval(ghost_map[i]);
    // Clear halo data
    ghost_map.clear();
    // Clear the ghost particle counter.
    _number_ghost = 0;
  }

  void SimData::removeHaloAndGhostParticles() {
    removeGhostParticles();
    removeHaloParticles();
  }

  int SimData::size() const {
    return _size;
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
    needs_remake = r;
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

  void SimData::quick_sort_help(int start, int end, int dim) {
    if (start<end && dim<sim_dimensions) {
      int partition = quick_sort_partition(start, end, dim);
      quick_sort_help(start, partition, dim);
      quick_sort_help(partition+1, end, dim);
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
        quick_sort_help(i*ds, (i+1)*ds, dim);
        recursion_help (i*ds, (i+1)*ds, dim+1);
      }
      // Potentially, there is some left over
      quick_sort_help(sort_bins*ds, _number-1, dim);
      recursion_help (sort_bins*ds, _number-1, dim+1);
    }
  }

}
