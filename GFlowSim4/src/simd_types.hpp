#ifndef __SIMD_TYPES_HPP__GFLOW__
#define __SIMD_TYPES_HPP__GFLOW__

#include "simd.hpp"

#if   SIMD_TYPE==SIMD_NONE
// The number of floats per vector
  const unsigned int simd_data_size = 1u;
  // What type we use as the vector of floats
  typedef float simd_float;

  // Functions
  inline simd_float simd_load (const float const *a)                   { return *a; }
  inline simd_float simd_set1 (const float a)                          { return a; }
  inline void       simd_store(const simd_float src, float *dest)      { *dest = src; }
  inline simd_float simd_add  (const simd_float a, const simd_float b) { return a+b; }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return a*b; }

#elif SIMD_TYPE==SIMD_SSE3
  // The number of floats per vector
  const unsigned int simd_data_size = 4u;
  // What type we use as the vector of floats
  typedef __m128 simd_float;

  // Functions
  inline simd_float simd_load (const float const *a)                   { return _mm_load_ps(a); }
  inline simd_float simd_set1 (const float a)                          { return _mm_set1_ps(a); }
  inline void       simd_store(const simd_float src, float *dest)      { _mm_store_ps(dest, src); }
  inline simd_float simd_add  (const simd_float a, const simd_float b) { return _mm_add_ps(a, b); }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return _mm_mul_ps(a, b); }

#elif SIMD_TYPE==SIMD_AVX or SIMD_TYPE==SIMD_AVX2
  // The number of floats per vector
  const unsigned int simd_data_size = 8u;
  // What type we use as the vector of floats
  typedef __m256 simd_float;

  // Functions
  inline simd_float simd_load (const float const *a)                   { return _mm256_load_ps(a); }
  inline simd_float simd_set1 (const float a)                          { return _mm256_set1_ps(a); }
  inline void       simd_store(const simd_float src, float *dest)      { _mm256_store_ps(dest, src); }
  inline simd_float simd_add  (const simd_float a, const simd_float b) { return _mm256_add_ps(a, b); }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return _mm256_mul_ps(a, b); }


#elif SIMD_TYPE==SIMD_MIC
  // The number of floats per vector
  const unsigned int simd_data_size = 16u;
  // What type we use as the vector of floats
  typedef __m512d simd_float;

  // Functions
  inline simd_float simd_load (const float const *a)                   { return _mm512_load_ps(a); }
  inline simd_float simd_set1 (const float a)                          { return _mm512_set1_ps(a); }
  inline void       simd_store(const simd_float src, float *dest)      { _mm512_store_ps(dest, src); }
  inline simd_float simd_add  (const simd_float a, const simd_float b) { _mm512_add_ps(dest, src); }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return _mm512_mul_ps(a, b); }

#endif

#if SIMD_TYPE!=SIMD_NONE
// Operators go here
#endif


#endif // __SIMD_TYPES_HPP__GFLOW__