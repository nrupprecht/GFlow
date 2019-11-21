// \file Template implementation file for particle-data-aos.hpp.

template<int dims> 
ParticleContainer_AOS<dims>::ParticleContainer_AOS() {
  // Essential vector data.
  vector_data_names = vector<string> {"X", "V", "F"};
  // Essential scalar data.
  scalar_data_names = vector<string> {"R", "Im"};
  // Essential integer data.
  integer_data_names = vector<string> {"Type", "ID"};
}

template<int dims> 
void ParticleContainer_AOS<dims>::initialize() {
  initialized = true;
}

template<int dims> 
bool ParticleContainer_AOS<dims>::is_soa() { 
  return false; 
}

template<int dims> 
void ParticleContainer_AOS<dims>::update() {

};

template<int dims> 
int ParticleContainer_AOS<dims>::add_vector_entry(const string& name) {
  if (!initialized && std::find(begin(vector_data_names), end(vector_data_names), name)==vector_data_names.end()) {
    vector_data_names.push_back(name);
    ++n_vectors;
    return n_vectors-1;
  }
  return -1;
}

template<int dims> 
int ParticleContainer_AOS<dims>::add_scalar_entry(const string& name) {
  if (!initialized && std::find(begin(scalar_data_names), end(scalar_data_names), name)==scalar_data_names.end()) {
    scalar_data_names.push_back(name);
    ++n_scalars;
    return n_scalars-1;
  }
  return -1;
}

template<int dims> 
int ParticleContainer_AOS<dims>::add_integer_entry(const string& name) {
  if (!initialized && std::find(begin(integer_data_names), end(integer_data_names), name)==integer_data_names.end()) {
    integer_data_names.push_back(name);
    ++n_integers;
    return n_integers-1;
  }
  return -1;
}

template<int dims> 
void ParticleContainer_AOS<dims>::clear_vec(int ar) {
  int address = ar*dims;
  for (int i=0; i<_size_owned; ++i, ar += data_width) 
    set1_vec<dims>(&data_ptr[address], 0);
}

template<int dims> 
int ParticleContainer_AOS<dims>::add_particle(vec<dims> x, real r, real mass, int type) {
  return add_particle(x, vec<dims>(), r, mass, type);
}

template<int dims> 
int ParticleContainer_AOS<dims>::add_particle(vec<dims> x, vec<dims> v, real r, real mass, int type) {
  if (_size_owned >= _capacity) {
    int new_capacity = _capacity + std::max(32, static_cast<int>(0.15*_capacity));
    resize(new_capacity);
  }

  // Add in particle by incrementing size counters, copying/clearing data.
  zeroVec(&data_ptr[data_width*_size_owned], data_width);
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
void ParticleContainer_AOS<dims>::reserve(uint s) {
  if (_capacity<s) resize(s);
}

template<int dims> 
int ParticleContainer_AOS<dims>::size() const {
  return _size_owned;
}

template<int dims> 
int ParticleContainer_AOS<dims>::number() const {
  return _number_owned;
}

template<int dims> 
void ParticleContainer_AOS<dims>::mark_for_removal(int id) {
  Type(id) = -1;
  remove_list.push_back(id);
}

template<int dims> 
void ParticleContainer_AOS<dims>::do_particle_removal() {
  
}

template<int dims> 
void ParticleContainer_AOS<dims>::resize(const int total_size) {
  if (total_size<=0) return;
  // In case we are resize to a smaller size.
  const int copy_size = std::min(_size_owned, total_size);
  // Create new pointer, copy data.
  real *new_data = new real[total_size];
  if (data_ptr) {
    delete [] data_ptr;
    copyVec(data_ptr, new_data, copy_size*data_width);
  }
  // Set ptr to new data.
  data_ptr = new_data;
  // Reset counters.
  _first_ghost = total_size;
  _number_ghost = 0;
  _size_ghost = 0;
  _capacity = total_size;
}
