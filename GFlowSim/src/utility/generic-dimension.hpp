#ifndef __GENERIC_DIMENSION_HPP__GFLOW__
#define __GENERIC_DIMENSION_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  inline bool isWrap(const BCFlag bc) {
    return bc == BCFlag::WRAP;
  }

  //! \brief Template function applying harmonic corrections to a distance vector.
  template<int d> inline void harmonic_correction(const BCFlag *bcs, RealType *r, const RealType *widths) {
    if (isWrap(bcs[d-1])) {
      RealType dX = widths[d-1] - fabs(r[d-1]);
      if (dX<fabs(r[d-1])) r[d-1] = r[d-1]>0 ? -dX : dX;
    } 
    // Recursion.
    harmonic_correction<d-1>(bcs, r, widths);
  }
  template<> inline void harmonic_correction<0>(const BCFlag *bcs, RealType *r, const RealType *widths) {};


  //! \brief Template function for the dot product of two vectors.
  template<int d> inline RealType dot_vec(const RealType *x, const RealType *y) {
    return x[d-1]*y[d-1] + dot_vec<d-1>(x, y);
  }
  template<> inline RealType dot_vec<1>(const RealType *x, const RealType *y) { return x[0]*y[0]; }; 


  //! \brief Template function for subtracting vectors.
  template<int d> inline void subtract_vec(const RealType *x, const RealType *y, RealType *z) {
    z[d-1] = x[d-1] - y[d-1];
    // Recursion.
    subtract_vec<d-1>(x, y, z);
  }
  template<> inline void subtract_vec<1>(const RealType *x, const RealType *y, RealType *z) {
    z[0] = x[0] - y[0];
  }

  
  //! \brief Template function for a vector *= a scalar.
  template<int d> inline void scalar_mult_eq_vec(RealType *x, const RealType m) {
    x[d-1] *= m;
    // Recursion.
    scalar_mult_eq_vec<d-1>(x, m);
  }
  template<> inline void scalar_mult_eq_vec<1>(RealType *x, const RealType m) {
    x[0] *= m;
  }

  
  //! \brief Template function for a vector += scalar * vector.
  template<int d> inline void sum_eq_vec_scaled(RealType *x, const RealType *y, const RealType v) {
    x[d-1] += v*y[d-1];
    // Recursion.
    sum_eq_vec_scaled<d-1>(x, y, v);
  }
  template<> inline void sum_eq_vec_scaled<1>(RealType *x, const RealType *y, const RealType v) {
    x[0] += v*y[0];
  }


  //! \brief Template function for copying (equating) a vector.
  template<int d> inline void copy_vec(const RealType *src, RealType *dst) {
    dst[d-1] = src[d-1];
    copy_vec<d-1>(src, dst);
  }
  template<> inline void copy_vec<1>(const RealType *src, RealType *dst) {
    dst[0] = src[0];
  }


}
#endif // __GENERIC_DIMENSION_HPP__GFLOW__
