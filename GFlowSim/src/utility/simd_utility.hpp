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

inline std::string simd_vec_to_vec_string(const simd_float *a, const int size) {
  stringstream stream;
  for (int s=0; s<simd_data_size; ++s) {
    stream << "{";
    for (int i=0; i<size; ++i) {
      stream << simd_get(s, a[i]);
      if (i!=size-1) stream << ",";
    }
    stream << "}";
    if (s!=simd_data_size-1) stream << ",";
  }
  std::string str;
  stream >> str;
  return str;
}

inline std::ostream& operator<<(std::ostream& out, simd_float x) {
  out << simd_to_str(x);
  return out;
}

//! @brief Store a vector in all simd coordinates
inline void simd_broadcast_vector(float *v, simd_float *sv, const int sim_dimensions) {
  for (int d=0; d<sim_dimensions; ++d) sv[d] = simd_set1(v[d]);
}

//! @brief Compute dot products
//!
//! { a_x1, a_x2, ... }, { a_y1, a_y2, ... }, ... * { a_x1, a_x2, ... }, { a_y1, a_y2, ... }, ... 
//! --> { a_x1 * b_x1 + a_y1*b_y1 + ... , a_x2 * b_x2 + a_y2 * b_y2 + ... , ... }
//!
//! So a[0] is a simd_float of all the x components of the first set of vectors, a[1] is a simd_float of all 
//! the y components of the first set of vectors, etc.
inline simd_float simd_dot_product(const simd_float *a, const simd_float *b, const int sim_dimensions) {
  // Accumulator
  simd_float acc = simd_zero;
  // Add components
  for (int i=0; i<sim_dimensions; ++i) {
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

inline void update_vector_data_size(const int *start_ids, float **vec, const simd_float *container, const int size, const int buffer_size) {
  for (int i=0; i<size; ++i) {
    int id = start_ids[i];
    for (int d=0; d<buffer_size; ++d)
      vec[id][d] += simd_get(i, container[d]);
  }
}

inline void simd_consolidate_update(float *v, const simd_float *vec, const int buffer_size) {
  for (int i=0; i<buffer_size; ++i) {
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

inline void simd_vector_sub(const simd_float *X1, const simd_float *X2, simd_float *dX, const int sim_dimensions) {
  for (int d=0; d<sim_dimensions; ++d)
    dX[d] = simd_sub(X1[d], X2[d]);
}

inline void simd_vector_sqr(const simd_float *X, simd_float& dX2, const int sim_dimensions) {
  dX2 = simd_set1(0.);
  for (int d=0; d<sim_dimensions; ++d)
    simd_plus_eq(dX2, simd_mult(X[d], X[d]));
}

#include "bounds.hpp"
using GFlowSimulation::Bounds;
using GFlowSimulation::BCFlag;

// Get the correct (minimal) displacement vector pointing from y to x
inline void simd_get_displacement(const simd_float *x, const simd_float *y, simd_float *dis, const Bounds B, 
  const BCFlag *boundaryConditions, int sim_dimensions) 
{
  for (int d=0; d<sim_dimensions; ++d) {
    dis[d] = simd_sub(x[d], y[d]);
    if (boundaryConditions[d]==BCFlag::WRAP) {

      simd_float dx = simd_sub(simd_set1(B.max[d] - B.min[d]), simd_abs(dis[d]));
      simd_float mask1 = simd_less_than(dx, simd_abs(dis[d]));
      simd_float inv_mask1 = simd_less_than(simd_abs(dis[d]), dx);
      simd_float mask2 = simd_less_than(simd_zero, dis[d]);

      // if (dx<fabs(dis[d])) dis[d] = dis[d]>0 ? -dx : dx;
      dis[d] = simd_mask(dis[d], inv_mask1);
      dis[d] = simd_add(
        simd_mask(simd_sub(dx, simd_mask(simd_mult(simd_set1(2.),dx), mask2)), mask1),
        dis[d]
      );
    }      
  }
}

#endif // __SIMD_UTILITY_HPP__GFLOW__
