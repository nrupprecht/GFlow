#ifndef __GENERIC_DIMENSION_HPP__GFLOW__
#define __GENERIC_DIMENSION_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  inline bool isWrap(const BCFlag bc) {
    return bc == BCFlag::WRAP;
  }

  //! \brief Template function applying harmonic corrections to a distance vector.
  template<int d> inline void harmonic_correction(const BCFlag *bcs, volatile RealType *r, const RealType *widths) {
    if (isWrap(bcs[d-1])) {
      RealType dX = widths[d-1] - fabs(r[d-1]);
      if (dX<fabs(r[d-1])) r[d-1] = r[d-1]>0 ? -dX : dX;
    } 
    // Recursion.
    harmonic_correction<d-1>(bcs, r, widths);
  }
  template<> inline void harmonic_correction<0>(const BCFlag *bcs, volatile RealType *r, const RealType *widths) {};

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
  template<int d> inline RealType dot_vec(const volatile RealType *x, const volatile RealType *y) {
    return x[d-1]*y[d-1] + dot_vec<d-1>(x, y);
  }
  template<> inline RealType dot_vec<1>(const volatile RealType *x, const volatile RealType *y) { return x[0]*y[0]; }; 

  //! \brief Template function for the dot product of two vectors.
  template<int d> inline RealType dot_vec(const RealType *x, const RealType *y) {
    return x[d-1]*y[d-1] + dot_vec<d-1>(x, y);
  }
  template<> inline RealType dot_vec<1>(const RealType *x, const RealType *y) { return x[0]*y[0]; }; 



  //! \brief Template function for the dot product of two vectors.
  template<int d> inline int dot_vec(const volatile int *x, const volatile int *y) {
    return x[d-1]*y[d-1] + dot_vec<d-1>(x, y);
  }
  template<> inline int dot_vec<1>(const volatile int *x, const volatile int *y) { return x[0]*y[0]; }; 

  //! \brief Template function for the dot product of two vectors.
  template<int d> inline int dot_vec(const int *x, const int *y) {
    return x[d-1]*y[d-1] + dot_vec<d-1>(x, y);
  }
  template<> inline int dot_vec<1>(const int *x, const int *y) { return x[0]*y[0]; }; 



  //! \brief Template function for adding vectors.
  template<int d> inline void add_vec(const volatile RealType *x, const volatile RealType *y, volatile RealType *z) {
    z[d-1] = x[d-1] + y[d-1];
    // Recursion.
    add_vec<d-1>(x, y, z);
  }
  template<> inline void add_vec<1>(const volatile RealType *x, const volatile RealType *y, volatile RealType *z) {
    z[0] = x[0] + y[0];
  }

  //! \brief Template function for adding vectors.
  template<int d> inline void add_vec(const RealType *x, const RealType *y, RealType *z) {
    z[d-1] = x[d-1] + y[d-1];
    // Recursion.
    add_vec<d-1>(x, y, z);
  }
  template<> inline void add_vec<1>(const RealType *x, const RealType *y, RealType *z) {
    z[0] = x[0] + y[0];
  }



  template<int d> inline void add_vec(const volatile int *x, const volatile int *y, volatile int *z) {
    z[d-1] = x[d-1] + y[d-1];
    // Recursion.
    add_vec<d-1>(x, y, z);
  }
  template<> inline void add_vec<1>(const volatile int *x, const volatile int *y, volatile int *z) {
    z[0] = x[0] + y[0];
  }

  template<int d> inline void add_vec(const int *x, const int *y, int *z) {
    z[d-1] = x[d-1] + y[d-1];
    // Recursion.
    add_vec<d-1>(x, y, z);
  }
  template<> inline void add_vec<1>(const int *x, const int *y, int *z) {
    z[0] = x[0] + y[0];
  }



  //! \brief Template function for subtracting vectors.
  template<int d> inline void subtract_vec(const volatile RealType *x, const volatile RealType *y, volatile RealType *z) {
    z[d-1] = x[d-1] - y[d-1];
    // Recursion.
    subtract_vec<d-1>(x, y, z);
  }
  template<> inline void subtract_vec<1>(const volatile RealType *x, const volatile RealType *y, volatile RealType *z) {
    z[0] = x[0] - y[0];
  }

  //! \brief Template function for subtracting vectors.
  template<int d> inline void subtract_vec(const RealType *x, const RealType *y, RealType *z) {
    z[d-1] = x[d-1] - y[d-1];
    // Recursion.
    subtract_vec<d-1>(x, y, z);
  }
  template<> inline void subtract_vec<1>(const RealType *x, const RealType *y, RealType *z) {
    z[0] = x[0] - y[0];
  }


  //! \brief Negate a vector.
  template<int d> inline void negate_vec(RealType *x) {
    x[d-1] = -x[d-1];
    negate_vec<d-1>(x);
  }
  template<> inline void negate_vec<1>(RealType *x) {
    x[0] = -x[0];
  };


  
  //! \brief Template function for a vector *= a scalar.
  template<int d> inline void scalar_mult_eq_vec(volatile RealType *x, const volatile RealType m) {
    x[d-1] *= m;
    // Recursion.
    scalar_mult_eq_vec<d-1>(x, m);
  }
  template<> inline void scalar_mult_eq_vec<1>(volatile RealType *x, const volatile RealType m) {
    x[0] *= m;
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



  template<int d> inline bool normalize_vec(RealType *x) {
    RealType rsqr = dot_vec<d>(x, x);
    if (rsqr==0) return false;
    scalar_mult_eq_vec<d>(x, 1./sqrt(rsqr));
    return true;
  }



  //! \brief Template function for a vector += vector.
  template<int d> inline void sum_eq_vec(volatile RealType *x, const volatile RealType *y) {
    x[d-1] += y[d-1];
    // Recursion.
    sum_eq_vec<d-1>(x, y);
  }
  template<> inline void sum_eq_vec<1>(volatile RealType *x, const volatile RealType *y) {
    x[0] += y[0];
  }

  //! \brief Template function for a vector += vector.
  template<int d> inline void sum_eq_vec(RealType *x, const RealType *y) {
    x[d-1] += y[d-1];
    // Recursion.
    sum_eq_vec<d-1>(x, y);
  }
  template<> inline void sum_eq_vec<1>(RealType *x, const RealType *y) {
    x[0] += y[0];
  }



  //! \brief Template function for a vector += vector.
  template<int d> inline void subtract_eq_vec(volatile RealType *x, const volatile RealType *y) {
    x[d-1] -= y[d-1];
    // Recursion.
    subtract_eq_vec<d-1>(x, y);
  }
  template<> inline void subtract_eq_vec<1>(volatile RealType *x, const volatile RealType *y) {
    x[0] -= y[0];
  }

  //! \brief Template function for a vector += vector.
  template<int d> inline void subtract_eq_vec(RealType *x, const RealType *y) {
    x[d-1] -= y[d-1];
    // Recursion.
    subtract_eq_vec<d-1>(x, y);
  }
  template<> inline void subtract_eq_vec<1>(RealType *x, const RealType *y) {
    x[0] -= y[0];
  }


  
  //! \brief Template function for a vector += scalar * vector.
  template<int d> inline void sum_eq_vec_scaled(volatile RealType *x, const volatile RealType *y, const RealType v) {
    x[d-1] += v*y[d-1];
    // Recursion.
    sum_eq_vec_scaled<d-1>(x, y, v);
  }
  template<> inline void sum_eq_vec_scaled<1>(volatile RealType *x, const volatile RealType *y, const RealType v) {
    x[0] += v*y[0];
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
  template<int d> inline void copy_vec(const volatile RealType *src, volatile RealType *dst) {
    dst[d-1] = src[d-1];
    copy_vec<d-1>(src, dst);
  }
  template<> inline void copy_vec<1>(const volatile RealType *src, volatile RealType *dst) {
    dst[0] = src[0];
  }

  //! \brief Template function for copying (equating) a vector.
  template<int d> inline void copy_vec(const RealType *src, RealType *dst) {
    dst[d-1] = src[d-1];
    copy_vec<d-1>(src, dst);
  }
  template<> inline void copy_vec<1>(const RealType *src, RealType *dst) {
    dst[0] = src[0];
  }



  template<int d> inline void copy_vec(const volatile int *src, volatile int *dst) {
    dst[d-1] = src[d-1];
    copy_vec<d-1>(src, dst);
  }
  template<> inline void copy_vec<1>(const volatile int *src, volatile int *dst) {
    dst[0] = src[0];
  }

  template<int d> inline void copy_vec(const int *src, int *dst) {
    dst[d-1] = src[d-1];
    copy_vec<d-1>(src, dst);
  }
  template<> inline void copy_vec<1>(const int *src, int *dst) {
    dst[0] = src[0];
  }



  //! \brief Set all the elements of a vector to be a single value
  template<int d> inline void set1_vec(int volatile *vec, int val) {
    vec[d-1] = val;
    set1_vec<d-1>(vec, val);
  }
  template<> inline void set1_vec<1>(int volatile *vec, int val) {
    vec[0] = val;
  }

  //! \brief Set all the elements of a vector to be a single value
  template<int d> inline void set1_vec(int *vec, int val) {
    vec[d-1] = val;
    set1_vec<d-1>(vec, val);
  }
  template<> inline void set1_vec<1>(int *vec, int val) {
    vec[0] = val;
  }



  //! \brief Set all the elements of a vector to be a single value
  template<int d> inline void set1_vec(RealType volatile *vec, RealType val) {
    vec[d-1] = val;
    set1_vec<d-1>(vec, val);
  }
  template<> inline void set1_vec<1>(RealType volatile *vec, RealType val) {
    vec[0] = val;
  }

  //! \brief Set all the elements of a vector to be a single value
  template<int d> inline void set1_vec(RealType *vec, RealType val) {
    vec[d-1] = val;
    set1_vec<d-1>(vec, val);
  }
  template<> inline void set1_vec<1>(RealType *vec, RealType val) {
    vec[0] = val;
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

  //! \brief Get the product of all the elements in an array.
  template<int d> inline int product(vector<int> array) {
    return array[d-1]*product<d-1>(array);
  }
  template<> inline int product<1>(vector<int> array) {
    return array[0];
  }



  //! Hadamard-equals operation
  template<int d> inline void hadamard_equals_vec(volatile RealType *x, const volatile RealType *y) {
    x[d-1] *= y[d-1];
    hadamard_equals_vec<d-1>(x, y);
  }
  template<> inline void hadamard_equals_vec<1>(volatile RealType *x, const volatile RealType *y) {
    x[0] *= y[0];
  }

  //! Hadamard-equals operation
  template<int d> inline void hadamard_equals_vec(RealType *x, const RealType *y) {
    x[d-1] *= y[d-1];
    hadamard_equals_vec<d-1>(x, y);
  }
  template<> inline void hadamard_equals_vec<1>(RealType *x, const RealType *y) {
    x[0] *= y[0];
  }



  // Helper functions for arbitrary for loop go into unnamed namespace.
  namespace {
    template<int d> inline void for_loop_helper(std::function<void(int*)> body, int *index, int *lower, int *upper) {
      for (index[d-1] = lower[d-1]; index[d-1]<upper[d-1]; ++index[d-1]) {
        for_loop_helper<d-1>(body, index, lower, upper);
      }
    }
    template<> inline void for_loop_helper<1>(std::function<void(int*)> body, int *index, int *lower, int *upper) {
      for (index[0] = lower[0]; index[0]<upper[0]; ++index[0]) {
        body(index);
      }
    }

    template<int d> inline void for_loop_helper_max(std::function<void(int*, int)> body, int *index, int *lower, int *upper, int& counts, int max_counts) {
      for (index[d-1] = lower[d-1]; index[d-1]<upper[d-1] && counts<max_counts; ++index[d-1]) {
        for_loop_helper_max<d-1>(body, index, lower, upper, counts, max_counts);
      }
    }
    template<> inline void for_loop_helper_max<1>(std::function<void(int*, int)> body, int *index, int *lower, int *upper, int& counts, int max_counts) {
      for (index[0] = lower[0]; index[0]<upper[0] && counts<max_counts; ++index[0], ++counts) {
        body(index, counts);
      }
    }
  }

  // Arbitrary dimensional for loop
  template<int d> inline void for_loop(std::function<void(int*)> body, int *lower, int *upper) {
    int index[d];
    copy_vec<d>(lower, index);
    // Use the for_loop_helper function with the index.
    for_loop_helper<d>(body, index, lower, upper);
  }

  // Arbitrary dimensional for loop that only executes max_counts steps.
  template<int d> inline void for_loop(std::function<void(int*, int)> body, int *lower, int *upper, int max_counts) {
    int index[d];
    copy_vec<d>(lower, index);
    int counts = 0;
    // Use the for_loop_helper function with the index.
    for_loop_helper_max<d>(body, index, lower, upper, counts, max_counts);
  }


  template<int d> constexpr inline int power(int x) {
    return x*power<d-1>(x);
  }
  template<> constexpr inline int power<1>(int x) {
    return x;
  }


}
#endif // __GENERIC_DIMENSION_HPP__GFLOW__
