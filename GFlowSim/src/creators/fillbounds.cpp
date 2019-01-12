#include "fillbounds.hpp"

namespace GFlowSimulation {

  RectangularBounds::RectangularBounds(int dim) : FillBounds(dim) {
    min = new RealType[dim];
    max = new RealType[dim];
  }

  RectangularBounds::RectangularBounds(Bounds bnds, int dim) : FillBounds(dim) {
    min = new RealType[dim];
    max = new RealType[dim];
    for (int d=0; d<dim; ++d) {
      min[d] = bnds.min[d];
      max[d] = bnds.max[d];
    }
  }

  RectangularBounds::~RectangularBounds() {
    if (min) delete [] min;
    if (max) delete [] max;
  }

  double RectangularBounds::vol() {
    float v = 1;
    for (int d=0; d<bnd_dimensions; ++d)
      v *= (max[d] - min[d]);
    return v;
  }

  void RectangularBounds::pick_position(RealType *x) {
    for (int d=0; d<bnd_dimensions; ++d) 
      x[d] = drand48()*(max[d] - min[d]) + min[d];
  }

  Bounds RectangularBounds::getBounds() {
    Bounds bnds(bnd_dimensions);
    for (int d=0; d<bnd_dimensions; ++d) {
      bnds.min[d] = min[d];
      bnds.max[d] = max[d];
    }
    return bnds;
  }

  SphericalBounds::SphericalBounds(int dim) : FillBounds(dim), radius(0) {
    center = new RealType[dim];
  }

  SphericalBounds::~SphericalBounds() {
    delete [] center;
  }

  double SphericalBounds::vol() {
    return sphere_volume(radius, bnd_dimensions);
  }

  void SphericalBounds::pick_position(RealType *x) {
    if (radius==0) {
      zeroVec(x, bnd_dimensions);
      return;
    }
    // Do this the dumb way for now, so it works in arbitrary # of dimensions
    bool good = false;
    while (!good) {
      for (int d=0; d<bnd_dimensions; ++d)
        x[d] = 2*(drand48() - 0.5)*radius + center[d];
      // Check whether the point is good
      good = sqr(x, bnd_dimensions)<=sqr(radius);
    }
  }

  Bounds SphericalBounds::getBounds() {
    Bounds bnds(bnd_dimensions);
    for (int d=0; d<bnd_dimensions; ++d) {
      bnds.min[d] = center[d] - radius;
      bnds.max[d] = center[d] + radius;
    }
    return bnds;
  }

}