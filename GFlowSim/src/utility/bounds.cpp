#include "bounds.hpp"
// Other files
#include "vectormath.hpp"

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

  //! @param m Bounds minima.
  //! @param M Bounds maxima.
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

  bool Bounds::operator==(const Bounds& b) const {
    if (dimensions!=b.dimensions) return false;
    for (int d=0; d<dimensions; ++d)
      if (min[d]!=b.min[d] || max[d]!=b.max[d]) return false;
    // The bounds passed all checks.
    return true;
  }

  bool Bounds::operator!=(const Bounds& b) const {
    return !(*this==b);
  }

  std::ostream& operator<<(std::ostream &out, Bounds bnds) {
    out << "{";
    for (int d=0; d<bnds.dimensions; ++d) {
      out << "{" << bnds.min[d] << ","  << bnds.max[d] << "}";
      if (d!=bnds.dimensions-1) out << ",";
    }
    out << "}";
    return out;
  }

  bool Bounds::contains(const RealType *x) const {
    // For the unlikely case where x[d] is exactly on a boundary,
    // we use half open [ ) boundaries.
    for (int d=0; d<dimensions; ++d)
      if (x[d]<min[d] || max[d]<=x[d]) return false;
    // All clear.
    return true;
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

  RealType Bounds::boundary() const {
    RealType v = vol(), bd = 0;
    for (int d=0; d<dimensions; ++d) {
      // Calculate the (perimeter, area, volume, etc) of each "face."
      // Two "faces," one on each side.
      bd += 2*v/(max[d] - min[d]);
    }
    return bd;
  }

  RealType Bounds::distance(const RealType *x) const {
    // Tally up distances in each dimension.
    RealType dsqr = 0;
    for (int d=0; d<dimensions; ++d) {
      if (x[d]<min[d]) dsqr += sqr(min[d]-x[d]);
      else if (max[d]<x[d]) dsqr += sqr(x[d]-max[d]);
    }
    // Return the distance.
    return sqrt(dsqr);
  }

  int Bounds::dims() const { 
    return dimensions; 
  }

  void Bounds::center(RealType *v) const {
    for (int d=0; d<dimensions; ++d) {
      v[d] = 0.5*(max[d] + min[d]);
    }
  }

  void Bounds::randomPoint(RealType *x) const {
    // We use drand because its fast and simple
    //! \todo Replace with a different RNG?
    for (int d=0; d<dimensions; ++d)
      x[d] = drand48()*wd(d) + min[d];
  }

  Bounds Bounds::intersection(const Bounds A, const Bounds B) {
    if (A.dimensions!=B.dimensions) throw BadDimension("Bounds intersection.");
    Bounds C(A.dimensions);
    bool null_intersection = false;
    for (int d=0; d<C.dimensions; ++d) {
      C.min[d] = GFlowSimulation::max(A.min[d], B.min[d]);
      C.max[d] = GFlowSimulation::min(A.max[d], B.max[d]);
      if (C.max[d]<=C.min[d]) {
        null_intersection = true;
        break;
      }
    }
    // If empty intersection, set to zeros.
    if (null_intersection) {
      for (int d=0; d<C.dimensions; ++d)
        C.min[d] = C.max[d] = 0.;
    }
    // Return the bounds.
    return C;
  }

  RealType max_width(const Bounds& bnds) {
    RealType a = bnds.wd(0), b = bnds.wd(1), c = bnds.wd(2);
    return max(a, b, c);
  }

  RealType min_width(const Bounds& bnds) {
    RealType a = bnds.wd(0), b = bnds.wd(1), c = bnds.wd(2);
    return min(a, b, c);
  }

}
