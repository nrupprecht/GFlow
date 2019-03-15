#include "particle-container.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  ParticleContainer::ParticleContainer(GFlow *gflow) : Base(gflow) {
    // Add default data entries
    addVectorData("X");
    addVectorData("V");
    addIntegerData("Type");
    addIntegerData("ID");
  };

  ParticleContainer::~ParticleContainer() {
    if (particle_data) delete [] particle_data;
    particle_data = nullptr;
  }

  //! \brief Initialize the object.
  void ParticleContainer::initialize() {
    // Call base's initialize function.
    Base::initialize();
    // The the offsets for the necessary data.
    _pos_offset = requestVectorData("X");
    _vel_offset = requestVectorData("V");
    _type_offset = requestIntegerData("Type");
    _id_offset = requestIntegerData("ID");
    // Offsets for unnecessary, but common, data
    _force_offset = requestVectorData("F");
    _sg_offset = requestScalarData("Sg");
    _im_offset = requestScalarData("Im");

    // Get ntypes from force master
    _ntypes = forceMaster->getNTypes();

    // Set flag
    _is_initialized = true;
  }

  void ParticleContainer::post_integrate() {
    // STUB
  }

  void ParticleContainer::reserve(int n_particles) {
    if (_size + n_particles >= _capacity) resize( n_particles - (_capacity - _size) );
  }

  void ParticleContainer::addParticle() {
    // Make sure there is enough space for the particle
    if (_size + 1 >= _capacity) resize(32);
    // Add type - default is type 0
    particle_data[_size*data_width + _type_offset] = byte_cast<float>(0);
    // Add id
    particle_data[_size*data_width + _type_offset] = byte_cast<float>(_next_global_id);

    // Increment number, id counter
    ++_number;
    ++_size;
    ++_next_global_id;
  }

  void ParticleContainer::addParticle(const RealType *x, const RealType *v, const RealType sg, const RealType im, const int type) {
    // Make sure there is enough space for the particle
    if (_size + 1 >= _capacity) resize(32);

    // Add position
    copyVec(x, &particle_data[_size*data_width + _pos_offset], 3);
    // Add velocity
    copyVec(v, &particle_data[_size*data_width + _vel_offset], 3);
    // Add type
    particle_data[_size*data_width + _type_offset] = byte_cast<float>(type);
    // Add id
    particle_data[_size*data_width + _type_offset] = byte_cast<float>(_next_global_id);

    // Add cutoff radius, if applicable
    if (_sg_offset>=0) particle_data[_size*data_width + _sg_offset] = sg;
    // Add cutoff inverse mass, if applicable
    if (_im_offset>=0) particle_data[_size*data_width + _im_offset] = im;

    // Increment number, id counter
    ++_number;
    ++_size;
    ++_next_global_id;
  }

  void ParticleContainer::markForRemoval(const int id) {
    if (Type(id)<0) return;
    // Mark for removal, clear some data
    remove_list.insert(id);
    setType(id, -1);
    /*
    id_map.erase(Id(id));
    Id(id) = -1;
    */
  }

  void ParticleContainer::doParticleRemoval() {

  }

  void ParticleContainer::exchangeParticles() {
    
  }

  void ParticleContainer::sortParticles() {
    
  }

  void ParticleContainer::updateHaloParticles() {

  }

  RealType* ParticleContainer::X(int i) {
    return &particle_data[i*data_width + _pos_offset];
  }

  RealType& ParticleContainer::X(int i, int d) {
    return particle_data[i*data_width + _pos_offset + d];
  }

  RealType* ParticleContainer::V(int i) {
    return &particle_data[i*data_width + _vel_offset];
  }

  RealType& ParticleContainer::V(int i, int d) {
    return particle_data[i*data_width + _vel_offset + d];
  }

  RealType* ParticleContainer::F(int i) {
    return &particle_data[i*data_width + _force_offset];
  }

  RealType& ParticleContainer::F(int i, int d) {
    return particle_data[i*data_width + _force_offset + d];
  }

  RealType& ParticleContainer::Sg(int i) {
    return particle_data[i*data_width + _sg_offset];
  }

  RealType& ParticleContainer::Im(int i) {
    return particle_data[i*data_width + _im_offset];
  }

  const RealType* ParticleContainer::X(int i) const {
    return &particle_data[i*data_width + _pos_offset];
  }

  const RealType* ParticleContainer::V(int i) const {
    return &particle_data[i*data_width + _vel_offset];
  }

  const RealType* ParticleContainer::F(int i) const {
    return &particle_data[i*data_width + _force_offset];
  }

  RealType ParticleContainer::Sg(int i) const {
    return particle_data[i*data_width + _sg_offset];
  }

  RealType ParticleContainer::Im(int i) const {
    return particle_data[i*data_width + _im_offset];
  }
  
  //! \brief Get the type of a particle. Since this information is encoded as a float, we cannot just return a reference to it.
  int ParticleContainer::Type(int i) {
    return byte_cast<int>(particle_data[i*data_width + _type_offset]);
  }

  void ParticleContainer::setType(int i, int type) {
    particle_data[i*data_width + _type_offset] = byte_cast<float>(type);
  }

  int ParticleContainer::Id(int i) {
    return byte_cast<int>(particle_data[i*data_width + _id_offset]);
  }

  RealType* ParticleContainer::particle(int i) {
    return &particle_data[i];
  }

  int ParticleContainer::size() {
    return _size;
  }

  int ParticleContainer::number() {
    return _number;
  }

  void ParticleContainer::addVectorData(const string& name) {
    // Only add an entry if it does not already exist.
    if (requestVectorData(name)==-1) {
      vector_data_offsets.insert(SIPair(name, vector_data_offsets.size()));
      data_width += sim_dimensions;
      _is_initialized = false;
    }
  }

  void ParticleContainer::addScalarData(const string& name) {
    // Only add an entry if it does not already exist.
    if (requestScalarData(name)==-1) {
      scalar_data_offsets.insert(SIPair(name, scalar_data_offsets.size()));
      ++data_width;
      _is_initialized = false;
    }
  }

  void ParticleContainer::addIntegerData(const string& name) {
    // Only add an entry if it does not already exist.
    if (requestIntegerData(name)==-1) {
      integer_data_offsets.insert(SIPair(name, integer_data_offsets.size()));
      ++data_width;
      _is_initialized = false;
    }
  }

  int ParticleContainer::requestVectorData(const string& name) {
    // Check if the data already exists
    auto it = vector_data_offsets.find(name);
    if (it!=vector_data_offsets.end()) return it->second;
    else return -1;
  }

  int ParticleContainer::requestScalarData(const string& name) {
    // Check if the data already exists
    auto it = scalar_data_offsets.find(name);
    if (it!=scalar_data_offsets.end()) return it->second;
    else return -1;
  }

  int ParticleContainer::requestIntegerData(const string& name) {
    // Check if the data already exists
    auto it = integer_data_offsets.find(name);
    if (it!=integer_data_offsets.end()) return it->second;
    else return -1;
  }

  void ParticleContainer::clearV() {
    // Clear all force vectors
    for (int i=0; i<_size; ++i)
      zeroVec(V(i), sim_dimensions);
  }

  void ParticleContainer::clearF() {
    // If this system does not have forces, return.
    if (_force_offset<0) return;
    // Clear all force vectors
    for (int i=0; i<_size; ++i) 
      zeroVec(F(i), sim_dimensions);
  }

  inline void ParticleContainer::resize(int additional_capacity) {
    // Only increase the size
    if (additional_capacity>0) {
      RealType *array = new RealType[(_capacity + additional_capacity)*data_width];
      copyVec(particle_data, array, _size);
      // Set type of remaining (new) particles slotes to be -1
      for (int i=_size; i<_capacity + additional_capacity; ++i) setType(i, -1);
      delete [] particle_data;
      particle_data = array;
      // Update capacity
      _capacity += additional_capacity;
    }
  }

  inline void ParticleContainer::swap_particle(int id1, int id2) {
    // Swap all entries of the particles.
    for (int i=0; i<data_width; ++i) 
      std::swap(particle(id1)[i], particle(id2)[i]);
  }

}