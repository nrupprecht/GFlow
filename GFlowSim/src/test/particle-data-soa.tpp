// \file Template implementation file for particle-data-soa.hpp.

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
int ParticleContainer_SOA<dims>::size() const {
  return _size_owned;
}

template<int dims>
int ParticleContainer_SOA<dims>::number() const {
  return _number_owned;
}

template<int dims>
void ParticleContainer_SOA<dims>::mark_for_removal(int id) {
  Type(id) = -1;
  remove_list.push_back(id);
}

template<int dims>
void do_particle_removal() {

}
