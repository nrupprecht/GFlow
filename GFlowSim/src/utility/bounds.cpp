#include "bounds.hpp"

namespace GFlowSimulation {

  Bounds::Bounds(const int dim) : dimensions(dim) {
    if (dim<=0) throw BadDimension();
    min = new RealType[dim]; 
    max = new RealType[dim];
    for (int d=0; d<dim; ++d) {
      min[d] = 0;
      max[d] = 0;
    }
  }

  Bounds::Bounds(const RealType *m, const RealType *M, const int dim) : dimensions(dim) {
    if (dim<=0) throw BadDimension();
    min = new RealType[dim];
    max = new RealType[dim];
    for (int d=0; d<dim; ++d) {
      min[d] = m[d];
      max[d] = M[d];
    }
  }

  Bounds::Bounds(const Bounds& b) : dimensions(b.dimensions) {
    // Set up arrays
    min = new RealType[dimensions];  
    max = new RealType[dimensions];
    // Copy bounds
    for (int d=0; d<dimensions; ++d) {
      min[d] = b.min[d];
      max[d] = b.max[d];
    }
  }

  Bounds::~Bounds() {
    if (min) delete [] min;
    if (max) delete [] max;
    min = nullptr;
    max = nullptr;
  }

  Bounds& Bounds::operator=(const Bounds& b) {
    if (dimensions!=b.dimensions) {
      if (min) delete [] min;
      if (max) delete [] max;
      dimensions = b.dimensions;
      min = new RealType[dimensions];  
      max = new RealType[dimensions];
    }
    else dimensions = b.dimensions;
    // Copy bounds
    for (int d=0; d<dimensions; ++d) {
      min[d] = b.min[d];
      max[d] = b.max[d];
    }
    return *this;
  }

  RealType Bounds::wd(int d) const { 
    return max[d] - min[d]; 
  }

  RealType Bounds::vol() const {
    RealType v = 1.0;
    for (int d=0; d<dimensions; ++d)
      v *= (max[d] - min[d]);
    return v;
  }

  int Bounds::dims() const { 
    return dimensions; 
  }

  void Bounds::center(RealType *v) const {
    for (int d=0; d<dimensions; ++d) {
      v[d] = 0.5*(max[d] + min[d]);
    }
  }

}