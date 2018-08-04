#ifndef __VECTOR_MATH_HPP__GFLOW__
#define __VECTOR_MATH_HPP__GFLOW__

#include "bounds.hpp"
#include "printingutility.hpp"

namespace GFlowSimulation {

  // Get the correct (minimal) displacement vector pointing from y to x
  inline void getDisplacement(const RealType *x, const RealType *y, RealType *dis, const Bounds B, const BCFlag *boundaryConditions) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) {
      dis[d] = x[d] - y[d];
      if (boundaryConditions[d]==BCFlag::WRAP) {
        RealType dx = B.max[d] - B.min[d] - fabs(dis[d]);
        if (dx<fabs(dis[d])) dis[d] = dis[d]>0 ? -dx : dx;
      }      
    }
  }

  template<int D=DIMENSIONS> inline RealType getDistanceSqrNoWrap(const RealType *x, const RealType *y) {
    RealType distance = 0;
    for (int d=0; d<D; ++d) distance += (x[d] - y[d]);
    return distance;
  }

  template<> inline RealType getDistanceSqrNoWrap<1>(const RealType *x, const RealType *y) {
    return x[0]-y[0];
  }

  template<> inline RealType getDistanceSqrNoWrap<2>(const RealType *x, const RealType *y) {
    return (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]);
  }

  template<> inline RealType getDistanceSqrNoWrap<3>(const RealType *x, const RealType *y) {
    return (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]);
  }

  template<typename T> T inline dotVec(const T *x, const T *y) {
    T dt(0);
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) dt += x[d]*y[d];
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
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) x[d] = T(0);
  }

  template<typename T> inline void setVec(T *x, const T val) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) x[d] = val;
  }

  template<typename T> inline void addVec(const T *x, const T *y, T *z) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) z[d] = x[d]+y[d];
  }

  template<typename T> inline void subtractVec(const T *x, const T *y, T *z) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) z[d] = x[d]-y[d];
  }

  // Times
  template<typename T> inline void scalarMultVec(const RealType scalar, const T *x, T *z) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) z[d] = scalar*x[d];
  }

  // Times equals
  template<typename T> inline void scalarMultVec(const RealType scalar, T *z) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) z[d] *= scalar;
  }

  template<typename T> inline void plusEqVec(T *z, const T *x) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) z[d] += x[d];
  }

  template<typename T> inline void minusEqVec(T *z, const T *x) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) z[d] -= x[d];
  }

  template<typename T> inline void copyVec(const T *x, T *target) {
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) target[d] = x[d];
  }

  template<typename T> inline void copyVec(const T *x, T *target, const unsigned int size) {
    for (int d=0; d<size; ++d) target[d] = x[d];
  }

  template<typename T> inline T magnitudeVec(T *vec) {
    return sqrt(sqr(vec));
  }

  template<typename T> inline void normalVec(const T *x, T *norm) {
    T mag = 0;
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) mag += sqr(x[d]);
    T invMag = mag>0 ? 1./sqrt(mag) : 0;
    scalarMultVec(invMag, x, norm);
  }

  template<typename T> inline void normalizeVec(T *norm) {
    T mag = 0;
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) mag += sqr(norm[d]);
    T invMag = mag>0 ? 1./sqrt(mag) : 0;
    scalarMultVec(invMag, norm);
  }

  inline void randomNormalVec(RealType *x) {
    // Random normal components
    // Note: Loop cannot be unrolled
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