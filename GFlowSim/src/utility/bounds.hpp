#ifndef __BOUNDS_HPP__GFLOW__
#define __BOUNDS_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  /*
  *  @struct BoundsPack
  *
  *  For use when the dimensions are not known, and therefore can not be hard-coded in (templates
  *  require constant values).
  *
  */
   struct BoundsPack {
    BoundsPack() {};

    BoundsPack(int dim) {
      min = vector<double>(dim, 0);
      max = vector<double>(dim, 0);
    }

    BoundsPack(RealType *m, RealType *M, int d) {
      min = vector<double>(d, 0);
      max = vector<double>(d, 0);
      for (int i=0; i<d; ++i) {
        min[i] = m[i];
        max[i] = M[i];
      }
    }

    //! @brief Get the widths in various dimensions of the bounds.
    double wd(int d) { return max[d] - min[d]; }
    //! @brief Get the dimensionality of the BoundsPack.
    int dims() { return min.size(); }

    vector<double> min, max;
  };

  /*
  *  @struct BoundsBase
  *
  *  Primary template to BoundsBase
  *  Can be specialized to Linear, Rectangular, Rectangular prismic, etc bounds.
  *  Uses full template specification
  *  For more info, see <https://en.cppreference.com/w/cpp/language/template_specialization>,
  *  <https://en.cppreference.com/w/cpp/language/partial_specialization>
  *
  */
  template <int D> struct BoundsBase {
    BoundsBase() {
      for (int d=0; d<D; ++d) {
        min[d] = 0; max[d] = 0;
      }
    }

    // Return the width of the bounds
    RealType wd(int d) {
      return max[d] - min[d];
    }

    RealType vol() {
      RealType V = 1.;
      for (int d=0; d<D; ++d) V *= (max[d] - min[d]);
      return V;
    }

    BoundsPack pack_up() {
      return BoundsPack(min, max, D);
    }

    RealType min[D], max[D];
  };

  // Zero dimensional bounds - empty class
  template<> struct BoundsBase<0> {};

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

    // Return the width of the bounds
    RealType wd(int) {
      return max[0] - min[0];
    }

    RealType vol() {
      return (max[0] - min[0]);
    }

    BoundsPack pack_up() {
      return BoundsPack(min, max, 1);
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

    // Return the width of the bounds
    RealType wd(int d) {
      return max[d] - min[d];
    }

    RealType vol() {
      return (max[0] - min[0])*(max[1] - min[1]);
    }

    BoundsPack pack_up() {
      return BoundsPack(min, max, 2);
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

    // Return the width of the bounds
    RealType wd(int d) {
      return max[d] - min[d];
    }

    RealType vol() {
      return (max[0] - min[0])*(max[1] - min[1])*(max[2] - min[2]);
    }

    BoundsPack pack_up() {
      return BoundsPack(min, max, 3);
    }

    // Data
    RealType min[3], max[3];
  };

  // Which bounds base we will primarily use --- call this simply [Bounds]
  typedef BoundsBase<DIMENSIONS> Bounds;

}

#endif // __BOUNDS_HPP__GFLOW__
