#ifndef __VECTOR_MATH_HPP__GFLOW__
#define __VECTOR_MATH_HPP__GFLOW__

#include "bounds.hpp"
#include "printingutility.hpp"

namespace GFlowSimulation {

  // rx = rx − length*floor(rx/length + 0.5);
  inline void getDisplacement(const RealType *x, const RealType *y, RealType *dis, const Bounds &B, const BCFlag *boundaryConditions, int dimensions) {
    for (int d=0; d<dimensions; ++d) {
      dis[d] = x[d] - y[d];
      if (boundaryConditions[d]==BCFlag::WRAP) {
        RealType dx = B.max[d] - B.min[d] - fabs(dis[d]);
        if (dx<fabs(dis[d])) dis[d] = dis[d]>0 ? -dx : dx;
      }      
    }
  }

  template<typename T> T inline dotVec(const T *x, const T *y, int dimensions) {
    T dt(0);
    for (int d=0; d<dimensions; ++d) dt += x[d]*y[d];
    return dt;
  }

  //! \brief A dot product function that can handle a more diverse set of inputs and outputs.
  template<typename T, typename U, typename V> T inline dotVec(const U *x, const V *y, int dimensions) {
    T dt(0);
    for (int d=0; d<dimensions; ++d) dt += x[d]*y[d];
    return dt;
  }

  template<typename T> void inline hadamardVec(const T *x, const T *y, T *z, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] = x[d]*y[d];
  }

  //! \brief Cross product - assumes 3 dimensions.
  template<typename T> inline void crossVec(const T *x, const T *y, T *z) {
    z[0] = x[1]*y[2]  - y[1]*x[2];
    z[1] = -x[0]*y[2] + y[0]*x[2];
    z[2] = x[0]*y[1]  - y[0]*x[1];
  }

  //! \brief Cross product of two 2d vectors. Answer is just a number (vector in the z direction)
  template<typename T> inline T crossVec2(const T *x, const T *y) {
    return x[0]*y[1] - x[1]*y[0];
  }

  //! \brief Triple product of two dimensional vectors, X x ( Y x Z).
  template<typename T> inline void tripleProduct2(const T *x, const T *y, const T *z, T *out) {
    RealType yz = y[0]*z[1] - y[1]*z[0];
    out[0] = x[1]*yz;
    out[1] = -x[0]*yz;
  }

  // Generic squaring function
  template<typename T> T inline sqr(const T x) {
    return x*x;
  }

  template<typename T> T inline sqr(const T *x, int dimensions) {
    T val(0);
    for (int d=0; d<dimensions; ++d) val += sqr(x[d]);
    return val;
  }

  inline RealType sqr(const RealType *x, int dimensions) {
    return dotVec(x, x, dimensions);
  }

  inline int sqr(const int *x, int dimensions) {
    return dotVec(x, x, dimensions);
  }

  inline RealType getDistanceSqrNoWrap(const RealType *x, const RealType *y, int dimensions) {
    RealType distance = 0;
    for (int d=0; d<dimensions; ++d) distance += sqr(x[d] - y[d]);
    return distance;
  }

  template<typename T> inline void zeroVec(T *x, int dimensions) {
    for (int d=0; d<dimensions; ++d) x[d] = T(0);
  }

  template<typename T> inline void setVec(T *x, const T val, int dimensions) {
    std::fill(x, x+dimensions, val);
  }

  template<typename T> inline void swapVec(T *x, T *y, int dimensions) {
    for (int d=0; d<dimensions; ++d) std::swap(x[d], y[d]);
  }

  template<typename T> inline void addVec(const T *x, const T *y, T *z, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] = x[d]+y[d];
  }

  template<typename T> inline void subtractVec(const T *x, const T *y, T *z, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] = x[d]-y[d];
  }

  // Times
  template<typename S, typename T> inline void scalarMultVec(const S scalar, const T *x, T *z, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] = scalar*x[d];
  }

  // Times equals
  template<typename S, typename T> inline void scalarMultVec(const S scalar, T *z, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] *= scalar;
  }

  template<typename T> inline void plusEqVec(T *z, const T *x, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] += x[d];
  }

  template<typename T> inline void plusEqVecScaled(T *z, const T *x, const T s, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] += s*x[d];
  }

  template<typename T> inline void minusEqVec(T *z, const T *x, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] -= x[d];
  }

  template<typename T> inline void minusEqVecScaled(T *z, const T *x, const T s, int dimensions) {
    for (int d=0; d<dimensions; ++d) z[d] -= s*x[d];
  }

  template<typename T> inline void copyVec(const T *x, T *target, int dimensions) {
    std::copy(x, x+dimensions, target);
  }

  template<typename T> inline void setVec(T *v, const int start, const int end, const T val) {
    for (int i=start; i<end; ++i) v[i] = val;
  }

  template<typename T> inline T magnitudeVec(T *vec, int dimensions) {
    return sqrt(sqr(vec, dimensions));
  }

  template<typename T> inline void negateVec(T *vec, int dimensions) {
    for (int d=0; d<dimensions; ++d) vec[d] = -vec[d];
  }

  template<typename T> inline void normalVec(const T *x, T *norm, int dimensions) {
    T mag = 0;
    for (int d=0; d<dimensions; ++d) mag += sqr(x[d]);
    T invMag = mag>0 ? 1./sqrt(mag) : 0;
    scalarMultVec(invMag, x, norm, dimensions);
  }

  template<typename T> inline void normalizeVec(T *norm, int dimensions) {
    T mag = 0;
    for (int d=0; d<dimensions; ++d) mag += norm[d]*norm[d];
    T invMag = mag>0 ? 1./sqrt(mag) : 0;
    scalarMultVec(invMag, norm, dimensions);
  }

  inline void randomNormalVec(RealType *x, int dimensions) {
    // Random normal components
    // Note: Loop cannot be unrolled
    for (int d=0; d<dimensions; ++d)
      x[d] = randNormal();
    // Normalize
    normalizeVec(x, dimensions);
  }

}
#endif // __VECTOR_MATH_HPP__GFLOW__
