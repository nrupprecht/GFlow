// \file Template implementation file for particle-data-aos.hpp.

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
