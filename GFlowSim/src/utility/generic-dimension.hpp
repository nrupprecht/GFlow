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



  template<int d> inline int get_index(const RealType *x, const RealType *min, const RealType *invW, const int *shift, const int *dims, const int *products) {
    // Calculate index.
    int index = static_cast<int>((x[d-1] - min[d-1])*invW[d-1]) + shift[d-1];
    // Keep in bounds.
    index = (dims[d-1]<=index) ? dims[d-1]-1 : index;
    index = (index<0) ? 0 : index;
    // Recursive.
    return index*products[d] + get_index<d-1>(x, min, invW, shift, dims, products);
  }
  template<> inline int get_index<1>(const RealType *x, const RealType *min, const RealType *invW, const int *shift, const int *dims, const int *products) {
    // Calculate index
    int index = static_cast<int>((x[0] - min[0])*invW[0]) + shift[0];
    // Keep in bounds.
    index = (dims[0]<=index) ? dims[0]-1 : index;
    index = (index<0) ? 0 : index;
    // Return
    return index*products[1];
  }



  // --- Simd related

  /*

  //! \brief Scatter a single vector across a buffer.
  template<int i> inline void scatter_single0(int *ids, RealType **v, RealType **buffer) {
    buffer[0][i] = v[ids[i]][0];
  }

  template<int d, int i> inline void scatter_single(int *ids, RealType **v, RealType **buffer) {
    // Move the (d-1)-th component.
    buffer[d-1][i] = v[ids[i]][d-1];
    // Recursion.
    if (d==1) scatter_single0<i>(ids, v, buffer);
    else scatter_single<d-1, i>(ids, v, buffer);
  }
  


  template<int d> inline void load_simd_float(RealType **buffer, simd_float *dst) {
    // Load the (d-1)-th entry in the buffer.
    dst[d-1] = simd_load(buffer[d-1]);
    // Recursion.
    load_simd_float<d-1>(buffer, dst);
  }
  template<> inline void load_simd_float<1>(RealType **buffer, simd_float *dst) {
    dst[0] = simd_load(buffer[0]);
  }


  //! \brief Scatter array of structures (individual vectors) contained in the vector pointer v into an structure of arrays
  //!  (arrays of vector components) contained in buffer.
  template<int d> inline void scatter0(int *ids, RealType **v, RealType **buffer) {
    scatter_single<d, 0>(ids, buffer, v);
  }
  template<int d, int w> inline void scatter(int *ids, RealType **v, RealType **buffer) {
    // Scatter the w-th particle.
    scatter_single<d, w-1>(ids, buffer, v);
    // Recursion.
    if (d==1) scatter0<d>(ids, buffer, v);
    else scatter<d, w-1>(ids, buffer, v);
  }


  template<int d, int w> inline void stride_buffer(RealType *buffer, RealType *stride[d]) {
    stride[d-1] = &buffer[(d-1)*w];
    // Recursion.
    if (w==1) stride[0] = &buffer[0];
    else stride_buffer<d-1, w>(buffer, stride);
  }

  //! \brief Load vectors at different positions (indicated by the ids array) in v into a simd_float vector, dst.
  template<int d> inline void load_vector(int *ids, RealType **v, simd_float *dst) {
    // Create a buffer for AOS to SOA conversion.
    RealType buffer[d*simd_data_size];
    RealType *buff[d];

    stride_buffer(buffer, buff);

    // Pack the v-data into the buffer.
    scatter<simd_data_size>(ids, buff, d);
    // Load the buffer's entries into a simd_float.
    load_simd_float<d>(buff, dst);
  }
  */


}
#endif // __GENERIC_DIMENSION_HPP__GFLOW__
