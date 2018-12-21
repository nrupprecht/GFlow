#include "simdata.hpp"
// Other files
#include "../utility/memory.hpp"

#include "integrator.hpp"

namespace GFlowSimulation {

  SimData::SimData(GFlow *gflow) : Base(gflow), bounds(Bounds(2)) {
    // Initialize vdata array
    vdata = vector<RealType**>(3, nullptr);
    // Put names into map
    vector_data_map.insert(SIPair("X", 0));
    vector_data_map.insert(SIPair("V", 1)); 
    vector_data_map.insert(SIPair("F", 2));
    // Initialize sdata array
    sdata = vector<RealType*>(2, nullptr);
    // Put names into map
    scalar_data_map.insert(SIPair("Sg", 0));
    scalar_data_map.insert(SIPair("Im", 1));
    // Initialize idata array
    idata = vector<int*>(2, nullptr);
    // Put names into map
    integer_data_map.insert(SIPair("Type", 0));
    integer_data_map.insert(SIPair("ID", 1));
    // Set up bounds to have the propper dimensions
    bounds = Bounds(sim_dimensions);
  }

  SimData::~SimData() {
    for (auto &v : vdata) {
      if (v) dealloc_array_2d(v);
      v = nullptr;
    }
    for (auto &s : sdata) {
      if (s) delete [] s;
      s = nullptr;
    }
    for (auto &i : idata) {
      if (i) delete [] i;
      i = nullptr;
    }
  }

  //! @brief Initialize the atom container.
  void SimData::initialize() {
    // Call base's initialize
    Base::initialize();

    // For now
    bounds = gflow->getBounds();

    // Sort the particles by position
    sortParticles();
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
    resetParticle(_size);
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
      resetParticle(_size);
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
    resetParticle(_size);
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
        move_particle(count_back, id);
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
    if (removed>0) needs_remake = true;
  }

  void SimData::exchangeParticles() {
    // @todo Implement
  }

  void SimData::sortParticles() {
    // Make sure all particles are valid, and compressed
    doParticleRemoval();
    // Quick sort
    quick_sort_help(0, _number-1, 0);
    recursion_help (0, _number-1, 1);
  }

  void SimData::updateHaloParticles() {
    for (int i=0; i<halo_map.size(); i+=2) {
      int hid = halo_map[i];
      int pid = halo_map[i+1];
      // Update force
      plusEqVec(vdata[2][pid], vdata[2][hid], sim_dimensions);
      // Clear halo particle force record
      zeroVec(vdata[2][hid], sim_dimensions);
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

  int SimData::request_vector_data(string name) {
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

  int SimData::request_scalar_data(string name) {
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

  int SimData::request_integer_data(string name) {
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

  int SimData::get_vector_data(string name) {
    // Check if the data already exists
    auto it = vector_data_map.find(name);
    if (it!=vector_data_map.end()) return it->second;
    else return -1;
  }

  int SimData::get_scalar_data(string name) {
    // Check if the data already exists
    auto it = scalar_data_map.find(name);
    if (it!=scalar_data_map.end()) return it->second;
    else return -1;
  }

  int SimData::get_integer_data(string name) {
    // Check if the data already exists
    auto it = integer_data_map.find(name);
    if (it!=integer_data_map.end()) return it->second;
    else return -1;
  }

  int SimData::size() const {
    return _number;
  }

  int SimData::number() const {
    return _number;
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
    return bounds;
  }

  bool SimData::getNeedsRemake() {
    return needs_remake;
  }

  void SimData::setNeedsRemake(bool r) {
    needs_remake = r;
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

  void SimData::resetParticle(int id) {
    for (auto v : vdata) zeroVec(v[id], sim_dimensions);
    for (auto s : sdata) s[id] = 0.;
    for (auto i : idata) i[id] = 0;
  }

  void SimData::move_particle(int src, int dst) {
    // Get global IDs
    int gs = Id(src);
    int gd = Id(dst);
    // If we are overwriting a valid particle
    if (Type(dst)>-1) {
      --_number;
      id_map.erase(gd);
    }
    // Transfer data
    for (auto v : vdata) copyVec(v[src], v[dst], sim_dimensions);
    for (auto s : sdata) s[dst] = s[src];
    for (auto i : idata) i[dst] = i[src];

    // Remap global id
    auto it = id_map.find(gs);
    if (id_map.end()!=it)
      it->second = dst;
    else throw false; // \todo Make an exception for this case.

    // We have invalidated the local id
    needs_remake = true;
  }

  void SimData::swap_particle(int id1, int id2) {
    // Get global IDs
    int g1 = Id(id1);
    int g2 = Id(id2);
    // Transfer data
    for (auto v : vdata) swapVec(v[id1], v[id2], sim_dimensions);
    for (auto s : sdata) std::swap(s[id1], s[id2]);
    for (auto i : idata) std::swap(i[id1], i[id2]);
    // Set global ids
    auto it = id_map.find(g1);
    if (id_map.end()!=it) it->second = id1;
    it = id_map.find(g2);
    if (id_map.end()!=it) it->second = id2;
    // Set flag
    needs_remake = true;
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
