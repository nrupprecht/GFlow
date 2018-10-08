#ifndef __SIMD_UTILITY_HPP__GFLOW__
#define __SIMD_UTILITY_HPP__GFLOW__

#include "simd_types.hpp"

// --- Functions that only depend on the simd names --- //

#include <string>
#include <sstream>

inline std::string simd_to_str(const simd_float a) {
  float b[simd_data_size];
  simd_store_u(a, b);
  // Create a string stream
  std::stringstream stream;
  std::string str;
  for (int i=0; i<simd_data_size; ++i) {
    stream << b[i];
    if (i!=simd_data_size-1) stream << ",";
  }
  stream >> str;
  // Return the string
  return str;
}

inline std::string simd_to_str(const simd_int a) {
  stringstream stream;
  for (int i=0; i<simd_data_size; ++i) {
    stream << simd_get_int(i, a);
    if (i!=simd_data_size-1) stream << ",";
  }
  string str;
  stream >> str;
  return str;
}

inline std::string simd_vec_to_str(const simd_float *a, const int size) {
  string str = "{";
  for (int i=0; i<size; ++i) {
    str += simd_to_str(a[i]);
    if (i!=size-1) str += "},{";
  }
  str += "}";
  return str;
}

//! @brief Store a vector in all simd coordinates
template<int dimensions> inline void simd_broadcast_vector(float *v, simd_float *sv) {
  for (int d=0; d<dimensions; ++d) sv[d] = simd_set1(v[d]);
}

//! @brief Compute dot products
//!
//! { a_x1, a_x2, ... }, { a_y1, a_y2, ... }, ... * { a_x1, a_x2, ... }, { a_y1, a_y2, ... }, ... 
//! --> { a_x1 * b_x1 + a_y1*b_y1 + ... , a_x2 * b_x2 + a_y2 * b_y2 + ... , ... }
//!
//! So a[0] is a simd_float of all the x components of the first set of vectors, a[1] is a simd_float of all 
//! the y components of the first set of vectors, etc.
template<int dimensions> inline simd_float simd_dot_product(const simd_float *a, const simd_float *b) {
  // Accumulator
  simd_float acc = simd_set_zero();
  // Add components
  for (int i=0; i<dimensions; ++i) {
    simd_plus_eq(acc, simd_mult(a[i], b[i]) );
  }
  return acc;
}

inline simd_float simd_mask_length(int length) {
  if (length>=simd_data_size) return simd_cast_float(simd_set1_int(simd_valid));
  else {
    simd_int mask = simd_set1_int(0);
    for (int i=0; i<length; ++i) {
      simd_set(i, mask, simd_valid);
    }
    return simd_cast_float(mask);
  }
}

inline void load_scalar_data_simd(const int *start_ids, float *array, simd_float &container, int size, float *temp) {
  // --- Put data into the temp array
  int i=0;
  for (; i<size; ++i) {
    int id = start_ids[i];
    temp[i] = array[id];
  }
  for (; i<simd_data_size; ++i) {
    temp[i] = 0.;
  }
  // --- Simd load
  container = simd_load_u(temp);
}

inline void load_vector_data_simd(const int *start_ids, float **vec, simd_float *container, int size, int sim_dimensions, float **temp) {
  // --- Put data into the temp array
  int i=0;
  for (; i<size; ++i) {
    int id = start_ids[i];
    for (int d=0; d<sim_dimensions; ++d)
      temp[d][i] = vec[id][d];
  }
  for (; i<simd_data_size; ++i) {
    for (int d=0; d<sim_dimensions; ++d)
      temp[d][i] = 0.;
  }
  // --- Simd load
  for (int d=0; d<sim_dimensions; ++d)
    container[d] = simd_load_u(temp[d]);
}

inline void update_vector_data_size(const int *start_ids, float **vec, simd_float *container, int size, int sim_dimensions) {
  for (int i=0; i<size; ++i) {
    int id = start_ids[i];
    for (int d=0; d<sim_dimensions; ++d)
      vec[id][d] += reinterpret_cast<float*>(&container[d])[i];
  }
}

inline void simd_consolidate_update(float *v, const simd_float *vec, const int size) {
  for (int i=0; i<size; ++i) {
    for (int sd=0; sd<simd_data_size; ++sd)
      v[i] += simd_get(sd, vec[i]);
  }
}
  
inline void simd_vector_set1(const float *V, simd_float *vec, int sim_dimensions) {
  for (int d=0; d<sim_dimensions; ++d)
    vec[d] = simd_set1(V[d]);
}

inline void simd_scalar_mult_vec(const simd_float scalar, const simd_float *vec, simd_float *out, const int sim_dimensions) {
  for (int d=0; d<sim_dimensions; ++d)
    out[d] = simd_mult(scalar, vec[d]);
}

inline void simd_vector_sub(const simd_float *X1, const simd_float *X2, simd_float *dX, int sim_dimensions) {
  for (int d=0; d<sim_dimensions; ++d)
    dX[d] = simd_sub(X1[d], X2[d]);
}

inline void simd_vector_sqr(const simd_float *X, simd_float& dX2, int sim_dimensions) {
  dX2 = simd_set1(0.);
  for (int d=0; d<sim_dimensions; ++d)
    simd_plus_eq(dX2, simd_mult(X[d], X[d]));
}

inline std::ostream& operator<<(std::ostream& out, simd_float x) {
  out << simd_to_str(x);
  return out;
}

#endif // __SIMD_UTILITY_HPP__GFLOW__
