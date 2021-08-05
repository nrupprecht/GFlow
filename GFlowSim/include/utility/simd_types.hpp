#ifndef __SIMD_TYPES_HPP__GFLOW__
#define __SIMD_TYPES_HPP__GFLOW__

#include "simd.hpp"

//! @brief Bitwise all true.
const int simd_valid = 0xffffffff;

#if SIMD_TYPE==SIMD_NONE
// The number of floats per vector
  const int simd_data_size = 1u;
  // What type we use as the vector of floats
  typedef float simd_float;
  // What type we use as the vector of integers
  typedef int simd_int;

  // Logic
  inline simd_int simd_and(const simd_int a, const simd_int b) { return a & b; }

  // Casting
  inline simd_int   simd_cast_int(const simd_float a) { return *reinterpret_cast<const int*>(&a); }
  inline simd_float simd_cast_float(const simd_int a) { return *reinterpret_cast<const float*>(&a); }
  inline void move_to_int(const simd_float source, simd_int& dest) { dest = static_cast<int>(source); }

  // Masking 
  inline simd_float simd_mask(const simd_float a, const simd_float m) { return simd_cast_float(simd_and(simd_cast_int(a), simd_cast_int(m))); }

  // Setting, loading, storing
  inline simd_float simd_load (const float *a)                         { return *a; }
  inline simd_float simd_load_u(const float *a)                        { return *a; }
  inline simd_float simd_set1 (const float a)                          { return a; }
  inline simd_int   simd_set1_int(const int a)                         { return a; }
  inline simd_float simd_set_zero ()                                   { return 0; }
  inline void       simd_store(const simd_float src, float *dest)      { *dest = src; }
  inline void       simd_store(const simd_int src, int *dest)          { *dest = src; }
  inline void       simd_store_u(const simd_float src, float *dest)    { *dest = src; }
  inline void       simd_store_u(const simd_int src, int *dest)        { *dest = src; }
  // Arithmetic
  inline simd_float simd_add  (const simd_float a, const simd_float b) { return a+b; }
  inline simd_float simd_sub  (const simd_float a, const simd_float b) { return a-b; }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return a*b; }
  inline simd_float simd_div(const simd_float a, const simd_float b)   { return a/b; }
  inline void    simd_plus_eq (simd_float &a, const simd_float b)      { a += b; }
  inline simd_float simd_sqrt(const simd_float a)                      { return sqrt(a); }
  inline simd_float simd_inv_sqrt(const simd_float a)                  { return 1./sqrt(a); }
  inline simd_float simd_abs(const simd_float a)                       { return fabs(a); }
  inline simd_float simd_negate(const simd_float a)                    { return -a; }

  // Comparisons
  inline simd_int simd_less_than(const simd_float a, const simd_float b) { return a<b ? 1 : 0; }

  // Special load
  template<int dimensions> inline simd_float simd_load_constant(const float *a, const int i) {
    // Need to set in "reverse" order
    return a[i/dimensions];
  }

  inline simd_float simd_load_constant(const float *a, const int i, int dimensions) {
    // Need to set in "reverse" order
    return a[i/dimensions];
  }

#elif SIMD_TYPE==SIMD_SSE3
  // The number of floats per vector
  const int simd_data_size = 4u;
  // What type we use as the vector of floats
  typedef __m128 simd_float;
  // What type we use as the vector of integers
  typedef __m128i simd_int;

  // Constants
  const simd_float SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));

  // Logic
  inline simd_int simd_and(const simd_int a, const simd_int b) { return _mm_and_si128(a, b); }

  // Casting
  inline simd_int   simd_cast_int(const simd_float a) { return _mm_castps_si128(a); }
  inline simd_float simd_cast_float(const simd_int a) { return _mm_castsi128_ps(a); }

  inline void move_to_int(const simd_float source, simd_int& dest) {
    // This is done serially.
    for (int i=0; i<simd_data_size; ++i)
      reinterpret_cast<int32_t*>(&dest)[i] = static_cast<int>(source[i]);
  }

  // Masking 
  inline simd_float simd_mask(const simd_float a, const simd_float m) { return simd_cast_float(simd_and(simd_cast_int(a), simd_cast_int(m))); }

  // Setting, loading, storing
  inline simd_float simd_load (const float *a)                         { return _mm_load_ps(a); }
  inline simd_float simd_load_u(const float *a)                        { return _mm_loadu_ps(a); }
  inline simd_float simd_set1 (const float a)                          { return _mm_set1_ps(a); }
  inline simd_int   simd_set1_int(const int a)                         { return _mm_set1_epi32(a); }
  inline simd_float simd_set_zero ()                                   { return _mm_setzero_ps(); }
  inline void       simd_store(const simd_float src, float *dest)      { _mm_store_ps(dest, src); }
  inline void       simd_store_u(const simd_float src, float *dest)    { _mm_storeu_ps(dest, src); }

   // Arithmetic
  inline simd_float simd_add  (const simd_float a, const simd_float b) { return _mm_add_ps(a, b); }
  inline simd_float simd_sub  (const simd_float a, const simd_float b) { return _mm_sub_ps(a, b); }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return _mm_mul_ps(a, b); }
  inline simd_float simd_div  (const simd_float a, const simd_float b) { return _mm_div_ps(a, b); }
  inline void       simd_plus_eq (simd_float &a, const simd_float b)   { a = _mm_add_ps(a, b); }
  inline simd_float simd_sqrt(const simd_float a)                      { return _mm_sqrt_ps(a); }
  inline simd_float simd_inv_sqrt(const simd_float a)                  { return _mm_rsqrt_ps(a); }
  inline simd_float simd_abs(const simd_float a)                       { return _mm_andnot_ps(SIGNMASK, a); }
  inline simd_float simd_negate(const simd_float a)                    { return _mm_xor_ps(a, SIGNMASK); }

  // Comparisons
  inline simd_float simd_less_than(const simd_float a, const simd_float b) { return _mm_cmplt_ps(a, b); }
  inline simd_float simd_clamp(const simd_float a) { return simd_mask(a, simd_less_than(simd_set1(0), a)); }
  inline simd_float simd_un_clamp(const simd_float a) { return simd_mask(a, simd_less_than(a, simd_set1(0))); }

  // Special load
  inline simd_float simd_load_constant(const float *a, const int i, int dimensions) {
    // Need to set in "reverse" order
    return _mm_set_ps(a[(i+3)/dimensions], a[(i+2)/dimensions], a[(i+1)/dimensions], a[i/dimensions]);
  }

  template<int dimensions> inline simd_float simd_load_constant(const float *a, const int i) {
    // Need to set in "reverse" order
    return _mm_set_ps(a[(i+3)/dimensions], a[(i+2)/dimensions], a[(i+1)/dimensions], a[i/dimensions]);
  }

  template<> inline simd_float simd_load_constant<1>(const float *a, const int i) {
    return simd_load_u(&a[i]);
  }

  template<> inline simd_float simd_load_constant<2>(const float *a, const int i) {
    // Need to set in "reverse" order
    return _mm_set_ps(a[i/2+1], a[i/2+1], a[i/2], a[i/2]);
  }

  template<> inline simd_float simd_load_constant<4>(const float *a, const int i) {
    return simd_set1(a[i/4]);
  }

  inline void simd_update_masked(simd_float& a, const simd_float b, const simd_float mask) {
    simd_float a2 = _mm_andnot_ps(mask, a);
    simd_float updates = _mm_and_ps(b, mask);
    a = a2 + updates;
  }


#elif SIMD_TYPE==SIMD_AVX or SIMD_TYPE==SIMD_AVX2
  // The number of floats per vector
  const int simd_data_size = 8u;
  // What type we use as the vector of floats
  typedef __m256 simd_float;
  // What type we use as the vector of integers
  typedef __m256i simd_int;

  // Constants
  const simd_float SIGNMASK = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));

  // Logic
  inline simd_int simd_and(const simd_int a, const simd_int b) { return _mm256_and_si256(a, b); }

  // Casting
  inline simd_int   simd_cast_int(const simd_float a) { return _mm256_castps_si256(a); }
  inline simd_float simd_cast_float(const simd_int a) { return _mm256_castsi256_ps(a); }

  inline void move_to_int(const simd_float source, simd_int& dest) {
    // This is done serially.
    for (int i=0; i<simd_data_size; ++i)
      reinterpret_cast<int32_t*>(&dest)[i] = static_cast<int>(source[i]);
  }

  // Masking 
  inline simd_float simd_mask(const simd_float a, const simd_float m) { return simd_cast_float(simd_and(simd_cast_int(a), simd_cast_int(m))); }

  // Setting, loading, storing
  inline simd_float simd_load (const float *a)                         { return _mm256_load_ps(a); }
  inline simd_float simd_load_u(const float *a)                        { return _mm256_loadu_ps(a); }
  inline simd_float simd_set1 (const float a)                          { return _mm256_set1_ps(a); }
  inline simd_int   simd_set1_int(const int a)                         { return _mm256_set1_epi32(a); }
  inline simd_float simd_set_zero ()                                   { return _mm256_setzero_ps(); }
  inline void       simd_store(const simd_float src, float *dest)      { _mm256_store_ps(dest, src); }
  inline void       simd_store_u(const simd_float src, float *dest)    { _mm256_storeu_ps(dest, src); }
  // Arithmetic
  inline simd_float simd_add  (const simd_float a, const simd_float b) { return _mm256_add_ps(a, b); }
  inline simd_float simd_sub  (const simd_float a, const simd_float b) { return _mm256_sub_ps(a, b); }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return _mm256_mul_ps(a, b); }
  inline simd_float simd_div  (const simd_float a, const simd_float b) { return _mm256_div_ps(a, b); }
  inline void       simd_plus_eq (simd_float &a, const simd_float b)   { a = _mm256_add_ps(a, b); }
  inline simd_float simd_sqrt(const simd_float a)                      { return _mm256_sqrt_ps(a); }
  inline simd_float simd_inv_sqrt(const simd_float a)                  { return _mm256_rsqrt_ps(a); }
  inline simd_float simd_abs(const simd_float a)                       { return _mm256_andnot_ps(SIGNMASK, a); }
  inline simd_float simd_negate(const simd_float a)                    { return _mm256_xor_ps(a, SIGNMASK); }
  
  // Comparisons
  inline simd_float simd_less_than(const simd_float a, const simd_float b) { return _mm256_cmp_ps(a, b, 1); }
  inline simd_float simd_clamp(const simd_float a) { return simd_mask(a, simd_less_than(simd_set1(0.), a)); }
  inline simd_float simd_un_clamp(const simd_float a) { return simd_mask(a, simd_less_than(a, simd_set1(0.))); }

  // Special load
  inline simd_float simd_load_constant(const float *a, const int i, int dimensions) {
    // Need to set in "reverse" order
    return _mm256_set_ps(a[(i+7)/dimensions], a[(i+6)/dimensions], a[(i+5)/dimensions], a[(i+4)/dimensions], 
      a[(i+3)/dimensions], a[(i+2)/dimensions], a[(i+1)/dimensions], a[i/dimensions]);
  }

  template<int dimensions> inline simd_float simd_load_constant(const float *a, const int i) {
    // Need to set in "reverse" order
    return _mm256_set_ps(a[(i+7)/dimensions], a[(i+6)/dimensions], a[(i+5)/dimensions], a[(i+4)/dimensions], 
      a[(i+3)/dimensions], a[(i+2)/dimensions], a[(i+1)/dimensions], a[i/dimensions]);
  }

  template<> inline simd_float simd_load_constant<1>(const float *a, const int i) {
    // Need to set in "reverse" order
    return simd_load_u(&a[i]);
  }

  template<> inline simd_float simd_load_constant<2>(const float *a, const int i) {
    // Need to set in "reverse" order
    return _mm256_set_ps(a[i/2+3], a[i/2+3], a[i/2+2], a[i/2+2], 
      a[i/2+1], a[i/2+1], a[i/2], a[i/2]);
  }

  template<> inline simd_float simd_load_constant<4>(const float *a, const int i) {
    // Need to set in "reverse" order
    return _mm256_set_ps(a[i/2+1], a[i/2+1], a[i/2+1], a[i/2+1], a[i/2], a[i/2], a[i/2], a[i/2]);
  }

  template<> inline simd_float simd_load_constant<8>(const float *a, const int i) {
    // Need to set in "reverse" order
    return simd_set1(a[i/8]);
  }

#elif SIMD_TYPE==SIMD_MIC
  // The number of floats per vector
  const int simd_data_size = 16u;
  // What type we use as the vector of floats
  typedef __m512d simd_float;
  // What type we use as the vector of integers
  typedef __m512i simd_int;

  // Logic
  inline simd_int simd_and(const simd_int a, const simd_int b) { return _mm512_and_si128(a, b); }

  // Casting
  inline simd_int   simd_cast_int(const simd_float a) { return _mm512_castps_si128(a); }
  inline simd_float simd_cast_float(const simd_int a) { return _mm512_castsi128_ps(a); }

  inline void move_to_int(const simd_float source, simd_int& dest) {
    // This is done serially.
    for (int i=0; i<simd_data_size; ++i)
      reinterpret_cast<int32_t*>(&dest)[i] = static_cast<int>(source[i]);
  }

  // Masking 
  inline simd_float simd_mask(const simd_float a, const simd_float m) { return simd_cast_float(simd_and(simd_cast_int(a), simd_cast_int(m))); }

  // Setting, loading, storing
  inline simd_float simd_load (const float const *a)                   { return _mm512_load_ps(a); }
  inline simd_float simd_load_u(const float const *a)                  { return _mm512_loadu_ps(a); }
  inline simd_float simd_set1 (const float a)                          { return _mm512_set1_ps(a); }
  inline simd_float simd_set_zero ()                                   { return _mm512_setzero_ps(); }
  inline void       simd_store(const simd_float src, float *dest)      { _mm512_store_ps(dest, src); }
//  inline void       simd_store(const simd_int src, int *dest)          { _mm512_store_si512(dest, src); }
  inline void       simd_store_u(const simd_float src, float *dest)    { _mm512_storeu_ps(dest, src); }
//  inline void       simd_store_u(const simd_int src, int *dest)        { _mm512_storeu_si512(dest, src); }
  // Arithmetic
  inline simd_float simd_add  (const simd_float a, const simd_float b) { return _mm512_add_ps(a, b); }
  inline simd_float simd_sub  (const simd_float a, const simd_float b) { return _mm512_sub_ps(a, b); }
  inline simd_float simd_mult (const simd_float a, const simd_float b) { return _mm512_mul_ps(a, b); }
  inline void    simd_plus_eq (simd_float &a, const simd_float b)      { a = _mm512_add_ps(a, b); }
  // Comparisons
  // inline simd_int simd_less_than(const simd_float a, const simd_float b) { return _mm512_castpd_si512(_mm512_cmp_pd(a, b, _CMP_LT_OS)); }

  // Special load
  inline simd_float simd_load_constant(const float *a, const int i, int dimensions) {
    // Need to set in "reverse" order
    return _mm512_set_ps(
      a[(i+15)/dimensions], a[(i+14)/dimensions], a[(i+13)/dimensions], a[(i+12)/dimensions], 
      a[(i+11)/dimensions], a[(i+10)/dimensions], a[(i+9)/dimensions], a[(i+8)/dimensions],
      a[(i+7)/dimensions], a[(i+6)/dimensions], a[(i+5)/dimensions], a[(i+4)/dimensions], 
      a[(i+3)/dimensions], a[(i+2)/dimensions], a[(i+1)/dimensions], a[i/dimensions]
    );
  }

  template<int dimensions> inline simd_float simd_load_constant(const float *a, const int i) {
    // Need to set in "reverse" order
    return _mm512_set_ps(
      a[(i+15)/dimensions], a[(i+14)/dimensions], a[(i+13)/dimensions], a[(i+12)/dimensions], 
      a[(i+11)/dimensions], a[(i+10)/dimensions], a[(i+9)/dimensions], a[(i+8)/dimensions],
      a[(i+7)/dimensions], a[(i+6)/dimensions], a[(i+5)/dimensions], a[(i+4)/dimensions], 
      a[(i+3)/dimensions], a[(i+2)/dimensions], a[(i+1)/dimensions], a[i/dimensions]
    );
  }
#endif

#endif // __SIMD_TYPES_HPP__GFLOW__
