// \file Template implementation file for particle-data-soa.hpp.

template<int dims>
ParticleContainer_SOA<dims>::ParticleContainer_SOA(GFlow *gflow) : Base(gflow) {
  // Essential vector data.
  vector_data = vector<real**> {nullptr, nullptr, nullptr};
  vector_data_names = vector<string> {"X", "V", "F"};
  // Essential scalar data.
  scalar_data = vector<real*> {nullptr, nullptr};
  scalar_data_names = vector<string> {"R", "Im"};
  // Essential integer data.
  integer_data = vector<int*> {nullptr, nullptr};
  integer_data_names = vector<string> {"Type", "ID"};
}

template<int dims>
ParticleContainer_SOA<dims>::~ParticleContainer_SOA() {
  for (auto& ptr : vector_data) 
    if (ptr) dealloc_array_2d(ptr);
  for (auto& ptr : scalar_data) 
    if (ptr) delete [] ptr;
  for (auto& ptr : integer_data) 
    if (ptr) delete [] ptr;
  // Clear vectors, just in case.
  vector_data.clear();
  scalar_data.clear();
  integer_data.clear();
}

template<int dims> 
void ParticleContainer_SOA<dims>::initialize() {
  initialized = true;
}

template<int dims> 
void ParticleContainer_SOA<dims>::pre_integrate() {
  // Clear the timed object timer.
  clear_timer();
}

template<int dims> 
void ParticleContainer_SOA<dims>::post_integrate() {

}

template<int dims>
bool ParticleContainer_SOA<dims>::is_soa() { 
  return true; 
}

template<int dims>
void ParticleContainer_SOA<dims>::update() {};

template<int dims>
int ParticleContainer_SOA<dims>::add_vector_entry(const string& name) {
  if (!initialized && std::find(begin(vector_data_names), end(vector_data_names), name)==vector_data_names.end()) {
    vector_data.push_back(nullptr);
    vector_data_names.push_back(name);
    ++n_vectors;
    return vector_data.size()-1;
  }
  return -1;
}

template<int dims>
int ParticleContainer_SOA<dims>::add_scalar_entry(const string& name) {
  if (!initialized && std::find(begin(scalar_data_names), end(scalar_data_names), name)==scalar_data_names.end()) {
    scalar_data.push_back(nullptr);
    scalar_data_names.push_back(name);
    ++n_scalars;
    return scalar_data.size()-1;
  }
  return -1;
}

template<int dims>
int ParticleContainer_SOA<dims>::add_integer_entry(const string& name) {
  if (!initialized && std::find(begin(integer_data_names), end(integer_data_names), name)==integer_data_names.end()) {
    integer_data.push_back(nullptr);
    integer_data_names.push_back(name);
    ++n_integers;
    return integer_data.size()-1;
  }
  return -1;
}

template<int dims>
void ParticleContainer_SOA<dims>::clear_vec(int ar) {
  real *data = static_cast<real*>(vector_data[ar][0]);
  for (int i=0; i<_size_owned*dims; ++i) data[i] = 0.;
}

template<int dims>
int ParticleContainer_SOA<dims>::add_particle(vec<dims> x, real r, real mass, int type) {
  return add_particle(x, vec<dims>(), r, mass, type);
}

template<int dims>
int ParticleContainer_SOA<dims>::add_particle(vec<dims> x, vec<dims> v, real r, real mass, int type) {
  if (_size_owned >= _capacity) {
    int new_capacity = _capacity + std::max(32, static_cast<int>(0.15*_capacity));
    resize(new_capacity);
  }
  // Add in particle by incrementing size counters, copying/clearing data.
  for (auto ptr : vector_data) vec<dims>(ptr[_size_owned]).zero();
  for (auto ptr : scalar_data) ptr[_size_owned] = static_cast<real>(0);
  for (auto ptr : integer_data) ptr[_size_owned] = 0;
  // Now set the data.
  X(_size_owned) = x;
  V(_size_owned) = v;
  R(_size_owned) = r;
  Im(_size_owned) = 1./mass;
  Type(_size_owned) = type;
  // Increment counters.
  ++_number_owned;
  ++_size_owned;
  // Return the address of the added particle.
  return _size_owned-1;
}

template<int dims>
void ParticleContainer_SOA<dims>::reserve(uint s) {
  if (_capacity<s) resize(s);
}

template<int dims>
void ParticleContainer_SOA<dims>::mark_for_removal(int id) {
  Type(id) = -1;
  remove_list.push_back(id);
}

template<int dims>
void ParticleContainer_SOA<dims>::do_particle_removal() {
  /*
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
  _number_owned -= remove_list.size();
  // Array is compressed.
  _size_owned = _number;
  // Clear the remove list
  remove_list.clear();

  // We need to update. Removed could be zero if, e.g. only and all of the N particles were removed.
  // Arguably, you still might want to remake. There could be extra entries in the verlet list that 
  // it would be better to get rid of.
  if (removed>0) set_needs_remake();

  // Start simdata timer.
  stop_timer();
  */
}

template<int dims>
void ParticleContainer_SOA<dims>::resize(const int total_size) {
  if (total_size<=0) return;
  // In case we are resize to a smaller size.
  const int copy_size = std::min(_size_owned, total_size);
  // Allocate space for vector data, and copy old data.
  for (auto& ptr : vector_data) {
    real **old_ptr = ptr;
    ptr = alloc_array_2d<real>(total_size, dims);
    if (copy_size*dims>0) copyVec(old_ptr, ptr, copy_size*dims);
    if (old_ptr) dealloc_array_2d(old_ptr);
  }
  // Allocate space for scalar data, and copy old data.
  for (auto& ptr : scalar_data) {
    real *old_ptr = ptr;
    ptr = new real[total_size];
    if (copy_size*dims>0) copyVec(old_ptr, ptr, copy_size);
    if (old_ptr) delete [] old_ptr;
  }
  // Allocate space for scalar data, and copy old data.
  for (auto& ptr : integer_data) {
    int *old_ptr = ptr;
    ptr = new int[total_size];
    if (copy_size*dims>0) copyVec(old_ptr, ptr, copy_size);
    if (old_ptr) delete [] old_ptr;
  }

  // Reset counters.
  _first_ghost = total_size;
  _number_ghost = 0;
  _size_ghost = 0;
  _capacity = total_size;
}
