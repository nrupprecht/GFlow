#ifndef __SIMD_GENERIC_HPP__GFLOW__
#define __SIMD_GENERIC_HPP__GFLOW__

#include "printingutility.hpp"
#include "simd_utility.hpp"
#include "vectormath.hpp"

using namespace GFlowSimulation;

template<typename T>
inline T set1(float a);
// Template specification
template<> inline float set1<float>(float a) { return a; }
template<> inline simd_float set1<simd_float>(float a) { return simd_set1(a); }

inline float      add(float a, float b) { return a+b; }
inline simd_float add(simd_float a, simd_float b) { return simd_add(a, b); }

inline float      sub(float a, float b) { return a-b; }
inline simd_float sub(simd_float a, simd_float b) { return simd_sub(a, b); }

inline float      mult(float a, float b) { return a*b; }
inline simd_float mult(simd_float a, simd_float b) { return simd_mult(a, b); }

inline float mask_value(float a, float m) { return a*m; }
inline simd_float mask_value(simd_float a, simd_float m) { return simd_mask(a, m); }

inline float clamp(float a) { return a<0 ? 0 : a; }
inline simd_float clamp(simd_float a) { return simd_clamp(a); }

inline void scalar_mult_vec(float scalar, const float *vec, float *out, int) {
  GFlowSimulation::scalarMultVec(scalar, vec, out);
}

inline void scalar_mult_vec(simd_float scalar, const simd_float *vec, simd_float *out, int sim_dimensions) { 
  simd_scalar_mult_vec(scalar, vec, out, sim_dimensions); 
}

inline void copy_negative(const float *src, float *dest, int sim_dimensions) {
  for (int i=0; i<sim_dimensions; ++i) dest[i] = -src[i];
}

inline void copy_negative(const simd_float *src, simd_float *dest, int sim_dimensions) {
  for (int i=0; i<sim_dimensions; ++i) dest[i] = -src[i];
}

inline std::string to_str(const float a) {
  return toStr(a); 
}

inline std::string to_str(const simd_float a) {
  return simd_to_str(a);
}

inline std::string to_str_vec(const float *a, int) {
  return PrintingUtility::toStrVec(a);
}

inline std::string to_str_vec(const simd_float *a, int sim_dimensions) {
  return simd_vec_to_vec_string(a, sim_dimensions);
}

//******************************
// Operators
//******************************

inline simd_float operator+(const simd_float a, const simd_float b) {
  return simd_add(a, b);
}

inline simd_float operator-(const simd_float a, const simd_float b) {
  return simd_sub(a, b);
}

inline simd_float operator*(const simd_float a, const simd_float b) {
  return simd_mult(a, b);
}

inline simd_float operator/(const simd_float a, const simd_float b) {
  return simd_div(a, b);
}

#endif // __SIMD_GENERIC_HPP__GFLOW__