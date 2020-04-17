

//! \brief Constexpr function that can be used in static asserts to make sure that the particle type is valid.
//! That way, if I decide to add more types later, I only have to change it here to effect all the asserts.
inline constexpr bool check_particle_type(unsigned pt) {
  return pt < max_particle_types;
}

// --- Particle insertion / deletion

template<unsigned particle_type> 
int SimData::addParticle(bool create_global_id) {
  // Just call the more general add particle function. That function has the assert, so there is no need for one here.
  return addParticle<particle_type>(1, create_global_id);
}

//! \param num The number of particle slots to add.
template<unsigned particle_type> 
int SimData::addParticle(int num, bool create_global_id) {
  static_assert(check_particle_type(particle_type));
  // Can't add a negative number of particles.
  if (num<=0) return -1;
  int capacity = data_entries[particle_type].capacity();
  int& size = _size[particle_type];
  if (size+num > capacity) {
    int additional_capacity = max(32, static_cast<int>(0.1*(num+size-capacity)));
    data_entries[particle_type].resize(additional_capacity);
  }
  int address = size;
  for (int i=0; i<num; ++i) {
    // Reset all data
    reset_particle<particle_type>(size);
    // Set type, give a global id
    Type<particle_type>(size) = 0;
    if (particle_type==0 && create_global_id) {
      //if (use_id_map) id_map[particle_type].emplace(size, next_global_id);
      if (use_id_map) id_map[particle_type].emplace(next_global_id, size);
      Id<particle_type>(size) = next_global_id;
      next_global_id += d_global_id;
    }
    ++_number[particle_type];
    ++size;
  }
  // Return first particle id.
  return address;
}

template<unsigned particle_type> 
int SimData::addParticle(const real *x, const real *v, const real sg, const real im, const int type) {
  static_assert(check_particle_type(particle_type));
  // If not enough spots to add a new owned particle, create more.
  int capacity = data_entries[particle_type].capacity();
  int& size = _size[particle_type];
  if (size+1 > capacity) {
    int additional_capacity = max(32, static_cast<int>(0.1*size));
    data_entries[particle_type].resize(additional_capacity);
  }
  // Reset all data
  reset_particle<particle_type>(size);
  // Set data
  copyVec(x, X<particle_type>(size), sim_dimensions);
  copyVec(v, V<particle_type>(size), sim_dimensions);
  Sg<particle_type>(size) = sg;
  Im<particle_type>(size) = im;
  Type<particle_type>(size) = type;
  if (particle_type==0) {
    //if (use_id_map) id_map[particle_type].emplace(size, next_global_id);
    if (use_id_map) id_map[particle_type].emplace(next_global_id, size);
    Id<0>(size) = next_global_id++;
    next_global_id += d_global_id;
  }
  ++_number[particle_type];
  ++size;
  // Return particle (local) id.
  return size-1;
}

template<unsigned particle_type> 
void SimData::markForRemoval(const int id) {
  static_assert(check_particle_type(particle_type));
  // If the particle has already been marked for removal, or is beyond the end of the array, return.
  if (Type<particle_type>(id)<0 || id>=_size[particle_type]) return;
  // Mark for removal, clear some data
  if (particle_type==0) remove_list.emplace(id); // All ghost particles will be removed eventually.
  if (use_id_map) remove_global_id(particle_type, Id<particle_type>(id));
  Type<particle_type>(id) = -1;
  // Set position to be far away, so it will not be able to apply force on real particles.
  // This is mostly for the case where a particle that is a ghost particle on another processor gets erased.
  // The fact that this happens is not communicated, but if we move the particle really far away, it will be
  // as if the particle disappeared.
  X<particle_type>(id, 0) = -1000000.f;
  // Decrement number of particles.
  _number[particle_type] -= 1;
}

// --- Data access

template<unsigned particle_type> 
vec_access SimData::X() {
  // Since position data must exist, we don't need checks.
  return VectorData<particle_type, false>(0);
}

template<unsigned particle_type> 
real* SimData::X(const int id) {
  static_assert(check_particle_type(particle_type));
  return X<particle_type>()(id);
}

template<unsigned particle_type> 
real& SimData::X(const int id, const int dim) {
  static_assert(check_particle_type(particle_type));
  return X<particle_type>()(id, dim);
}

template<unsigned particle_type> 
vec_access SimData::V() {
  // Since velocity data must exist, we don't need checks.
  return VectorData<particle_type, false>(1);
}

template<unsigned particle_type> 
real* SimData::V(const int id) {
  static_assert(check_particle_type(particle_type));
  return V<particle_type>()(id);
}

template<unsigned particle_type> 
real& SimData::V(const int id, const int dim) {
  static_assert(check_particle_type(particle_type));
  return V<particle_type>()(id, dim);
}

template<unsigned particle_type> 
vec_access SimData::F() {
  // Since force data must exist, we don't need checks.
  return VectorData<particle_type, false>(2);
}

template<unsigned particle_type> 
real* SimData::F(const int id) {
  static_assert(check_particle_type(particle_type));
  return F<particle_type>()(id);
}

template<unsigned particle_type> 
real& SimData::F(const int id, const int dim) {
  static_assert(check_particle_type(particle_type));
  return F<particle_type>()(id, dim);
}

// This is the main function for getting vector data.
template<unsigned particle_type, bool safe_checks> 
vec_access SimData::VectorData(const int entry) {
  static_assert(check_particle_type(particle_type));
  // If we are using safe checks, make sure entry corresponds to a valid entry.
  if (safe_checks) {
    if (0<=entry && entry<nvectors()) return data_entries[particle_type].get_vdata(entry); 
    else return vec_access();
  }
  // If we are not using safe checks, assume we are safe.
  else return data_entries[particle_type].get_vdata(entry); 
}

template<unsigned particle_type> 
vec_access SimData::VectorData(const string& name) {
  const int entry = getVectorData(name);
  return VectorData<particle_type>(entry);
}

template<unsigned particle_type> 
real* SimData::VectorData(const int entry, const int id) {
  return VectorData<particle_type>(entry)(id);
}

// --- Scalar data ---

template<unsigned particle_type> 
scalar_access SimData::Sg() {
  return ScalarData<particle_type>(0);
}

template<unsigned particle_type> 
real& SimData::Sg(int id) {
  return Sg<particle_type>()(id);
}

template<unsigned particle_type> 
scalar_access SimData::Im() {
  return ScalarData<particle_type>(1);
}

template<unsigned particle_type> real& 
SimData::Im(int id) {
  return Im<particle_type>()(id);
}

// This is the main function for getting scalar data.
template<unsigned particle_type, bool safe_checks> 
scalar_access SimData::ScalarData(const int entry) {
  static_assert(check_particle_type(particle_type));
  // If we are using safe checks, make sure entry corresponds to a valid entry.
  if (safe_checks) {
    if (0<=entry && entry<nscalars()) return data_entries[particle_type].get_sdata(entry); 
    else return scalar_access();
  }
  // If we are not using safe checks, assume we are safe.
  else return data_entries[particle_type].get_sdata(entry); 
}

template<unsigned particle_type, bool safe_checks> 
scalar_access SimData::ScalarData(const string& name) {
  const int entry = getScalarData(name);
  return ScalarData<particle_type, safe_checks>(entry);
}

template<unsigned particle_type> 
real& SimData::ScalarData(const int entry, const int id) {
  return ScalarData<particle_type>(entry)(id);
}

// --- Integer data ---

template<unsigned particle_type> 
integer_access SimData::Type() {
  // Since type data must exist, we don't need checks.
  return IntegerData<particle_type, false>(0);
}

template<unsigned particle_type> 
int& SimData::Type(const int id) {
  return Type<particle_type>()(id);
}

template<unsigned particle_type> 
integer_access SimData::Id() {
// Since global id data must exist, we don't need checks.
  return IntegerData<particle_type, false>(1);
}

template<unsigned particle_type> 
int& SimData::Id(const int id) {
  return Id<particle_type>()(id);
}

// This is the main function for getting integer data.
template<unsigned particle_type, bool safe_checks> 
integer_access SimData::IntegerData(const int entry) {
  static_assert(check_particle_type(particle_type));
  // If we are using safe checks, make sure entry corresponds to a valid entry.
  if (safe_checks) {
    if (0<=entry && entry<nintegers()) return data_entries[particle_type].get_idata(entry); 
    else return integer_access();
  }
  // If we are not using safe checks, assume we are safe.
  else return data_entries[particle_type].get_idata(entry); 
}

template<unsigned particle_type> 
integer_access SimData::IntegerData(const string& name) {
  const int entry = getScalarData(name);
  return IntegerData<particle_type>(entry);
}

template<unsigned particle_type> 
int& SimData::IntegerData(const int entry, const int id) {
  return IntegerData<particle_type>(entry)(id);
}

template<unsigned particle_type>
void SimData::pack_buffer_relative(const vector<int>& id_list, vector<real>& buffer, const Vec& center_point) {
  int size = id_list.size();
  if (buffer.size()<size*data_width) buffer.resize(size*data_width);
  // Send the actual data. Copy data into buffer.
  RealType dx[4]; // Assumes sim_dimensions <= 4.
  int n_vectors = nvectors(), n_scalars = nscalars(), n_integers = nintegers();
  auto x = X<particle_type>();
  for (int j=0; j<size; ++j) {
    int id = id_list[j];
    // Get the position of the particle, relative to the other processor.
    gflow->getDisplacement(x(id), center_point.data, dx);
    plusEqVec(dx, center_point.data, sim_dimensions);
    // Copy particle information to the buffer, using the relative position. 
    // \todo Automate a way to specify arbitrary subsets of the particle data to send.
    copyVec(dx, &buffer[data_width*j], sim_dimensions); // Position
    // Send the rest of the data the normal way. Pack vector data.
    for (int i=1; i<n_vectors; ++i) 
      copyVec(VectorData<particle_type>(i, id), &buffer[data_width*j + i*sim_dimensions], sim_dimensions);
    // Pack scalar data.
    for (int i=0; i<n_scalars; ++i) 
      buffer[data_width*j + n_vectors*sim_dimensions + i] = ScalarData<particle_type>(i, id);
    // Pack integer data.
    for (int i=0; i<n_integers; ++i)
      buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i] = byte_cast<RealType>(IntegerData<particle_type>(i, id));
  }
  
}

template<unsigned particle_type>
void SimData::unpack_buffer(const int n_particles, const vector<real>& buffer) {
  int n_vectors = nvectors(), n_scalars = nscalars(), n_integers = nintegers();
  for (int j=0; j<n_particles; ++j) {
    // Add a spot for a particle, then copy the data into this particle.
    int id = addParticle<particle_type>(false); // Get a local id for new particle.
    // Unpack vector data.
    for (int i=0; i<n_vectors; ++i) 
      copyVec(&buffer[data_width*j + i*sim_dimensions], VectorData<particle_type>(i, id), sim_dimensions);
    // Unpack scalar data.
    for (int i=0; i<n_scalars; ++i) 
      ScalarData<particle_type>(i, id) = buffer[data_width*j + n_vectors*sim_dimensions + i];
    // Unpack integer data.
    for (int i=0; i<n_integers; ++i)
      IntegerData<particle_type>(i, id) = byte_cast<int>(buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i]);
  } 
}

// --- Helper functions

template<unsigned particle_type> 
void SimData::reset_particle(int id) {
  for (int i=0; i<nvectors(); ++i) zeroVec(VectorData<particle_type>(i, id), sim_dimensions);
  for (int i=0; i<nscalars(); ++i) ScalarData<particle_type>(i, id) = 0.;
  for (int i=0; i<nintegers(); ++i) IntegerData<particle_type>(i, id) = -1;
}
