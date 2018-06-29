#ifndef __VECTOR_MATH_HPP__GFLOW__
#define __VECTOR_MATH_HPP__GFLOW__

#include "bounds.hpp"
#include "printingutility.hpp"

namespace GFlowSimulation {

  // Get the correct (minimal) displacement vector pointing from y to x
 inline void getDisplacement(const RealType *x, const RealType *y, RealType *dis, const Bounds B, const bool *wrap) {
    for (int d=0; d<DIMENSIONS; ++d) {
      dis[d] = x[d] - y[d];
      if (wrap[d]) {
        RealType dx = B.max[d] - B.min[d] - fabs(dis[d]);
        if (dx<fabs(dis[d])) dis[d] = dis[d]>0 ? -dx : dx;
      }
    }
  }

  template<typename T> T inline dotVec(const T *x, const T *y) {
    T dt(0);
    for (int i=0; i<DIMENSIONS; ++i) dt += x[i]*y[i];
    return dt;
  }

  // Generic squaring function
  template<typename T> T inline sqr(const T x) {
    return x*x;
  }

  // Vector squaring function
  inline RealType sqr(RealType *x) {
    return dotVec(x,x);
  }

  template<typename T> inline void zeroVec(T *x) {
    for (int d=0; d<DIMENSIONS; ++d) x[d] = T(0);
  }

  template<typename T> inline void addVec(const T *x, const T *y, T *z) {
    for (int d=0; d<DIMENSIONS; ++d) z[d] = x[d]+y[d];
  }

  template<typename T> inline void subtractVec(const T *x, const T *y, T *z) {
    for (int d=0; d<DIMENSIONS; ++d) z[d] = x[d]-y[d];
  }

  // Times
  template<typename T> inline void scalarMultVec(const RealType scalar, const T *x, T *z) {
    for (int d=0; d<DIMENSIONS; ++d) z[d] = scalar*x[d];
  }

  // Times equals
  template<typename T> inline void scalarMultVec(const RealType scalar, T *z) {
    for (int d=0; d<DIMENSIONS; ++d) z[d] *= scalar;
  }

  template<typename T> inline void plusEqVec(T *z, const T *x) {
    for (int d=0; d<DIMENSIONS; ++d) z[d] += x[d];
  }

  template<typename T> inline void minusEqVec(T *z, const T *x) {
    for (int d=0; d<DIMENSIONS; ++d) z[d] -= x[d];
  }

  template<typename T> inline void copyVec(const T *x, T *target) {
    for (int d=0; d<DIMENSIONS; ++d) target[d] = x[d];
  }

  template<typename T> inline T magnitudeVec(T *vec) {
    return sqrt(sqr(vec));
  }

  template<typename T> inline void normalVec(const T *x, T *norm) {
    T mag = 0;
    for (int d=0; d<DIMENSIONS; ++d) mag += sqr(x[d]);
    T invMag = 1./sqrt(mag);
    scalarMultVec(invMag, x, norm);
  }

  template<typename T> inline void normalizeVec(T *norm) {
    T mag = 0;
    for (int d=0; d<DIMENSIONS; ++d) mag += sqr(norm[d]);
    T invMag = 1./sqrt(mag);
    scalarMultVec(invMag, norm);
  }

  inline void randomNormalVec(RealType *x) {
    // Random normal components
    for (int d=0; d<DIMENSIONS; ++d)
      x[d] = randNormal();
    // Normalize
    normalizeVec(x);
  }

  inline string toStrVec(RealType *vec) {
    return ("{"+PrintingUtility::toStrVec(vec)+"}");
  }

  inline string toStrVec(int *vec) {
    return ("{"+PrintingUtility::toStrVec(vec)+"}");
  }

}
#endif // __VECTOR_MATH_HPP__GFLOW__