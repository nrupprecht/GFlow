#ifndef __BOUNDS_HPP__GFLOW__
#define __BOUNDS_HPP__GFLOW__

namespace GFlowSimulation {

  struct Bounds {
    // Default constructor
    Bounds() : left(0), right(0), bottom(0), top(0) {};
    // Setting constructor
    Bounds(RealType l, RealType r, RealType b, RealType t) : left(l), right(r), bottom(b), top(t) {};

    // The data
    RealType left, right, bottom, top;
  };

}

#endif // __BOUNDS_HPP__GFLOW__