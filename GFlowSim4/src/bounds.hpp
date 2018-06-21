#ifndef __BOUNDS_HPP__GFLOW__
#define __BOUNDS_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  /*
  *  @struct BoundsBase
  *  Can be specialized to Linear, Rectangular, Rectangular prismic, etc bounds.
  *  Uses full template specification
  *  For more info, see <https://en.cppreference.com/w/cpp/language/template_specialization>
  *
  */
  template <int D> struct BoundsBase {};

  // One dimensional bounds base
  template<> struct BoundsBase<1> {
    // Default constructor
    BoundsBase() {
      min[0] = 0; max[0] = 0;
    }
    // Setting constructor
    BoundsBase(RealType left, RealType right) {
      min[0] = left; max[0] = right;
    }

    // Data
    RealType min[1], max[1];
  };

  // Two dimensional bounds base
  template<> struct BoundsBase<2> {
    // Default constructor
    BoundsBase() {
      min[0] = 0; max[0] = 0;
      min[1] = 0; max[1] = 0;
    }
    // Setting constructor
    BoundsBase(RealType left, RealType right, RealType bottom, RealType top) {
      min[0] = left;   max[0] = right;
      min[1] = bottom; max[1] = top;
    }

    // Data
    RealType min[2], max[2];
  };

  // Three dimensional bounds base
  template<> struct BoundsBase<3> {
    // Default constructor
    BoundsBase() {
      min[0] = 0; max[0] = 0;
      min[1] = 0; max[1] = 0;
      min[2] = 0; max[2] = 0;
    }
    // Setting constructor
    BoundsBase(RealType left, RealType right, RealType bottom, RealType top, RealType close, RealType far) {
      min[0] = left;   max[0] = right;
      min[1] = bottom; max[1] = top;
      min[2] = close;  max[2] = far;
    }

    // Data
    RealType min[3], max[3];
  };

  // Other
  // ...

  // Which bounds base we will use --- call this simply [Bounds]
  typedef BoundsBase<DIMENSIONS> Bounds;

}

#endif // __BOUNDS_HPP__GFLOW__