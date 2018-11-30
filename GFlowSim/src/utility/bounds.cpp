#include "bounds.hpp"

namespace GFlowSimulation {

  BoundsPack::BoundsPack() : dimensions(0), min(nullptr), max(nullptr) {};

  BoundsPack::BoundsPack(const int dim) : dimensions(dim) {
    min = new RealType[dim]; 
    max = new RealType[dim];
    for (int d=0; d<dim; ++d) {
      min[d] = 0;
      max[d] = 0;
    }
  }

  BoundsPack::BoundsPack(const RealType *m, const RealType *M, const int d) : dimensions(d) {
    min = new RealType[d];
    max = new RealType[d];
    for (int i=0; i<d; ++i) {
      min[i] = m[i];
      max[i] = M[i];
    }
  }

  BoundsPack::~BoundsPack() {
    if (min) delete [] min;
    if (max) delete [] max;
  }

  //! @brief Get the widths in various dimensions of the bounds.
  RealType BoundsPack::wd(int d) { return max[d] - min[d]; }

  //! @brief Get the dimensionality of the BoundsPack.
  int BoundsPack::dims() { return dimensions; }

  //! @brief Set the input vector to be the center of the bounds.
  void BoundsPack::center(RealType *v) {
    for (int i=0; i<dimensions; ++i) {
      v[i] = 0.5*(min[i] - max[i]);
    }
  }

}