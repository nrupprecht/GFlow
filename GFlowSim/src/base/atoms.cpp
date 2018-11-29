#include "atoms.hpp"
// Other files
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  Atoms::Atoms(GFlow *gflow) : Base(gflow) {
    // Initialize vdata array
    vdata = vector<RealType**>(3, nullptr);
    // Initialize sdata array
    sdata = vector<RealType*>(2, nullptr);
    // Initialize idata array
    idata = vector<int*>(2, nullptr);
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
    number = num;
    last_owned = num;
  }

  void Atoms::addParticle(const RealType *x, const RealType *v, const RealType sg, const RealType im, const int type) {
    // If not enough spots to add a new owned particle, create more
    if (number < last_owned) resize_owned(32);
    copyVec(x, X(number), sim_dimensions);
    copyVec(v, V(number), sim_dimensions);
    Sg(number) = sg;
    Im(number) = im;
    Type(number) = type;
    Id(number) = next_global_id++;
    ++number;
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

  void Atoms::resize_owned(int num) {
    

  }

}