#include "region.hpp"
// Other files
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  // --- Rectangular region

  RectangularRegion::RectangularRegion(int d) : Region(d), bounds(d) {};

  RectangularRegion::RectangularRegion(int d, const Bounds& b) : Region(d), bounds(b) {};

  void RectangularRegion::pick_position(RealType *x) {
    bounds.randomPoint(x);
  }

  bool RectangularRegion::contains(RealType *x) {
    return bounds.contains(x);
  }

  RealType RectangularRegion::vol() {
    return bounds.vol();
  }

  Bounds RectangularRegion::get_bounding_box() {
    return bounds;
  }

  RealType& RectangularRegion::min(int d) {
    return bounds.min[d];
  }

  RealType& RectangularRegion::max(int d ) {
    return bounds.max[d];
  }

  void RectangularRegion::setBounds(const Bounds& b) {
    bounds = b;
  }

  // --- Spherical region

  SphericalRegion::SphericalRegion(int d) : Region(d), _radius(0) {
    _center = new RealType[dimensions];
    zeroVec(_center, dimensions);
  }

  SphericalRegion::SphericalRegion(int d, const RealType *c, const RealType r) : Region(d), _radius(r) {
    _center = new RealType[dimensions];
    copyVec(c, _center, dimensions);
  }

  SphericalRegion::~SphericalRegion() {
    delete [] _center;
  }

  void SphericalRegion::pick_position(RealType *x) {
    if (_radius==0) {
      zeroVec(x, dimensions);
      return;
    }
    // Do this the naieve way for now, so it works in arbitrary # of dimensions
    bool good = false;
    Bounds bnds = get_bounding_box();
    while (!good) {
      bnds.randomPoint(x);
      // Check whether the point is good
      good = contains(x);
    }
  }

  bool SphericalRegion::contains(RealType *x) {
    RealType rsqr = 0;
    for (int d=0; d<dimensions; ++d) rsqr += sqr(x[d] - _center[d]);
    // Check whether the point is inside the sphere
    return rsqr<sqr(_radius);
  }

  RealType SphericalRegion::vol() {
    return sphere_volume(_radius, dimensions);
  }

  Bounds SphericalRegion::get_bounding_box() {
    // Create the smallest bounds that encloses the sphere
    Bounds bnds(dimensions);
    for (int d=0; d<dimensions; ++d) {
      bnds.min[d] = _center[d] - _radius;
      bnds.max[d] = _center[d] + _radius;
    }
    // Return the bounds
    return bnds;
  }

  RealType& SphericalRegion::center(int d) {
    return _center[d];
  }

  RealType& SphericalRegion::radius() {
    return _radius;
  }

  void SphericalRegion::setCenter(const RealType *r) {
    copyVec(r, _center, dimensions);
  }

  void SphericalRegion::setRadius(const RealType r) {
    _radius = r;
  }

}
