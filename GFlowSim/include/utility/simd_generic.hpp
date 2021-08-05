#ifndef __SIMD_GENERIC_HPP__GFLOW__
#define __SIMD_GENERIC_HPP__GFLOW__

#include "printingutility.hpp"
#include "simd_utility.hpp"
#include "vectormath.hpp"

using namespace GFlowSimulation;

template<typename T> inline T set1(float a);
// Template specification
template<> inline float set1<float>(const float a) { return a; }

inline float add(const float a, const float b) { return a+b; }

inline float sub(const float a, const float b) { return a-b; }

inline void sub(const float *a, const float *b, float *c, int sim_dimensions) {
  for (int d=0; d<sim_dimensions; ++d) c[d] = a[d] - b[d];
}

inline float mult(const float a, const float b) { return a*b; }

inline float dot(const float *a, const float *b, const int sim_dimensions) {
  return dotVec(a, b, sim_dimensions);
}

inline float mask_value(float a, float m) { return a*m; }

inline float clamp(float a) { return a<0 ? 0 : a; }

inline float un_clamp(float a) { return a<0 ? a : 0; }

inline void scalar_mult_vec(float scalar, const float *vec, float *out, int sim_dimensions) {
  GFlowSimulation::scalarMultVec(scalar, vec, out, sim_dimensions);
}

inline void copy_negative(const float *src, float *dest, int sim_dimensions) {
  for (int i=0; i<sim_dimensions; ++i) dest[i] = -src[i];
}

inline std::string to_str(const float a) {
  return toStr(a); 
}

inline std::string to_str_vec(const float *a, int sim_dimensions) {
  return toStrVec(a, sim_dimensions);
}

#if SIMD_TYPE!=SIMD_NONE
// Template specification
template<> inline simd_float set1<simd_float>(const float a) { return simd_set1(a); }
inline simd_float add(const simd_float a, const simd_float b) { return simd_add(a, b); }
inline simd_float sub(const simd_float a, const simd_float b) { return simd_sub(a, b); }
inline void       sub(const simd_float *a, const simd_float *b, simd_float *c, int sim_dimensions) {
  simd_vector_sub(a, b, c, sim_dimensions);
}
inline simd_float mult(simd_float a, simd_float b) { return simd_mult(a, b); }
inline simd_float dot(const simd_float *a, const simd_float *b, const int sim_dimensions) {
  return simd_dot_product(a, b, sim_dimensions);
}
inline simd_float mask_value(simd_float a, simd_float m) { return simd_mask(a, m); }
inline simd_float clamp(simd_float a) { return simd_clamp(a); }
inline simd_float un_clamp(simd_float a) { return simd_un_clamp(a); }
inline void scalar_mult_vec(simd_float scalar, const simd_float *vec, simd_float *out, int sim_dimensions) { 
  simd_scalar_mult_vec(scalar, vec, out, sim_dimensions); 
}
inline void copy_negative(const simd_float *src, simd_float *dest, int sim_dimensions) {
  for (int i=0; i<sim_dimensions; ++i) dest[i] = -src[i];
}
inline std::string to_str(const simd_float a) {
  return simd_to_str(a);
}
inline std::string to_str_vec(const simd_float *a, int sim_dimensions) {
  return simd_vec_to_vec_string(a, sim_dimensions);
}
#endif

#endif // __SIMD_GENERIC_HPP__GFLOW__
