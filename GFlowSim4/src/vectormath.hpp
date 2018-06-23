#ifndef __VECTOR_MATH_HPP__GFLOW__
#define __VECTOR_MATH_HPP__GFLOW__

#include "bounds.hpp"

namespace GFlowSimulation {

  // Get the correct (minimal) displacement vector pointing from y to x
 inline void getDisplacement(RealType *x, RealType *y, RealType *dis, const Bounds B, const bool *wrap) {
    for (int d=0; d<DIMENSIONS; ++d) {
      dis[d] = x[d] - y[d];
      if (wrap[d]) {
        RealType dx = B.max[d] - B.min[d] - fabs(dis[d]);
        if (dx<fabs(dis[d])) dis[d] = dis[d]>0 ? -dx : dx;
      }
    }
  }

  template<typename T> T dotVec(T *x, T *y) {
    T dt(0);
    for (int i=0; i<DIMENSIONS; ++i) dt += x[i]*y[i];
    return dt;
  }

  // Generic squaring function
  template<typename T> T sqr(T x) {
    return x*x;
  }

  // Vector squaring function
  RealType sqr(RealType *x) {
    return dotVec(x,x);
  }

  template<typename T> void zeroVec(T *x) {
    for (int d=0; d<DIMENSIONS; ++d)
      x[d] = T(0);
  }

  template<typename T> void addVec(T *x, T *y, T *z) {
    for (int d=0; d<DIMENSIONS; ++d) z[d] = x[d]+y[d];
  }

}
#endif // __VECTOR_MATH_HPP__GFLOW__