#include "atoms.hpp"
// Other files
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  Atoms::Atoms(GFlow *gflow) : Base(gflow), bounds(Bounds(2)) {
    // Initialize vdata array
    vdata = vector<RealType**>(3, nullptr);
    // Initialize sdata array
    sdata = vector<RealType*>(2, nullptr);
    // Initialize idata array
    idata = vector<int*>(2, nullptr);
    // Set up bounds to have the propper dimensions
    bounds = Bounds(sim_dimensions);
  }

  Atoms::~Atoms() {
    for (auto &v : vdata)
      dealloc_array_2d(v);
    for (auto &s : sdata)
      if (s) delete [] s;
    for (auto &i : idata)
      if (i) delete [] i;
  }

  //! @brief Initialize the atom container.
  void Atoms::initialize() {
    Base::initialize();

    // For now
    bounds = gflow->getBounds();
  }

  //! @brief Reserve space for particles, extending the lengths of all arrays to the requested size.
  void Atoms::reserve(int num) {
    for (auto &v : vdata) {
      if (v) dealloc_array_2d(v);
      alloc_array_2d<RealType>(num, sim_dimensions);
    }
    for (auto &s : sdata) {
      if (s) delete [] s;
      s = new RealType[num];
    }
    for (auto &i : idata) {
      if (i) delete [] i;
      i = new int[num];
    }
    // Reset numbers
    number = 0;
    last_valid = 0;
    owned_capacity = num;
    last_extended = num;
    capacity = num;
  }

  void Atoms::addParticle(const RealType *x, const RealType *v, const RealType sg, const RealType im, const int type) {
    // If not enough spots to add a new owned particle, create more
    if (number < owned_capacity) resize_owned(32);
    copyVec(x, X(number), sim_dimensions);
    copyVec(v, V(number), sim_dimensions);
    Sg(number) = sg;
    Im(number) = im;
    Type(number) = type;
    Id(number) = next_global_id++;
    ++number;
    ++last_valid;
  }

  void Atoms::markForRemoval(const int id) {
    // @todo Implement
  }

  void Atoms::doParticleRemoval() {
    // @todo Implement
  }

  void Atoms::exchangeParticles() {
    // @todo Implement
  }

  void Atoms::updateHaloParticles() {
    for (int i=0; i<halo_map.size(); i+=2) {
      int hid = halo_map[i];
      int pid = halo_map[i+1];
      // Update force
      plusEqVec(vdata[2][pid], vdata[2][hid], sim_dimensions);
      // Clear halo particle force record
      zeroVec(vdata[2][hid], sim_dimensions);
    }
  }

  RealType** Atoms::X() {
    return vdata[0];
  }

  RealType* Atoms::X_arr() {
    return *vdata[0];
  }

  RealType* Atoms::X(int i) {
    return vdata[0][i];
  }

  RealType& Atoms::X(int i, int d) {
    return vdata[0][i][d];
  }

  RealType** Atoms::V() {
    return vdata[1];
  }

  RealType* Atoms::V_arr() {
    return *vdata[0];
  }

  RealType* Atoms::V(int i) {
    return vdata[0][i];
  }

  RealType& Atoms::V(int i, int d) {
    return vdata[1][i][d];
  }

  RealType** Atoms::F() {
    return vdata[2];
  }

  RealType* Atoms::F_arr() {
    return *vdata[2];
  }

  RealType* Atoms::F(int i) {
    return vdata[2][i];
  }

  RealType& Atoms::F(int i, int d) {
    return vdata[2][i][d];
  }

  RealType* Atoms::Sg() {
    return sdata[0];
  }

  RealType& Atoms::Sg(int i) {
    return sdata[0][i];
  }

  RealType* Atoms::Im() {
    return sdata[1];
  }

  RealType& Atoms::Im(int i) {
    return sdata[1][i];
  }

  int* Atoms::Type() {
    return idata[0];
  }

  int& Atoms::Type(int i) {
    return idata[0][i];
  }

  int* Atoms::Id() {
    return idata[1];
  }

  int& Atoms::Id(int i) {
    return idata[1][i];
  }

  int Atoms::size() {
    return number;
  }

  int Atoms::Number() {
    return number;
  }

  int Atoms::getLocalID(int global) const {
    auto it = id_map.find(global);
    // Return the global iterator. We use -1 to mean "no such particle."
    return it==id_map.end() ? -1 : it->second;
  }

  int Atoms::getNextGlobalID() const {
    return next_global_id;
  }

  const BCFlag* Atoms::getBCs() const {
    return Base::gflow->getBCs();
  }

  Bounds Atoms::getBounds() const {
    return bounds;
  }

  void Atoms::resize_owned(int num) {
    // Compute new capacity
    int new_capacity = capacity + num;
    // Allocate new vector data arrays
    for (auto &v : vdata) {
      RealType **nv = alloc_array_2d<RealType>(new_capacity, sim_dimensions);
      // Transfer data
      copyData(v, nv, 0, last_valid, 0); // Copy owned data
      copyData(v, nv, last_valid+1, last_extended, last_valid+num); // Copy extended data
      // Delete old array, set new
      delete [] v;
      v = nv;
    }
    // Allocate new scalar data arrays
    for (auto &s : sdata) {
      RealType *ns = new RealType[new_capacity];
      // Transfer data
      copyData(s, ns, 0, last_valid, 0); // Copy owned data
      copyData(s, ns, last_valid+1, last_extended, last_valid+num); // Copy extended data
      // Delete old array, set new
      delete [] s;
      s = ns;
    }
    // Allodate new integer data
    // Allocate new scalar data arrays
    for (auto &i : idata) {
      int *ni = new int[new_capacity];
      // Transfer data
      copyData(i, ni, 0, last_valid, 0); // Copy owned data
      copyData(i, ni, last_valid+1, last_extended, last_valid+num); // Copy extended data
      // Delete old array, set new
      delete [] i;
      i = ni;
    }

    // Set new sizes
    owned_capacity += num;
    last_extended += num;
    capacity += num;;
  }

  void Atoms::copyData(RealType **source, RealType **dest, int start, int end, int new_start) {
    int j = 0;
    for (int i=start; i<end; ++i, ++j)
      copyVec(source[i], dest[new_start+j], sim_dimensions);
  }

  void Atoms::copyData(RealType *source, RealType *dest, int start, int end, int new_start) {
    int j = 0;
    for (int i=start; i<end; ++i, ++j)
      dest[new_start+j] = source[i];
  }

  void Atoms::copyData(int *source, int *dest, int start, int end, int new_start) {
    int j = 0;
    for (int i=start; i<end; ++i, ++j)
      dest[new_start+j] = source[i];
  }

}