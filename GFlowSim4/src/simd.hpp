#pragma once

#ifndef __SIMD_HPP__GFLOW__
#define __SIMD_HPP__GFLOW__


// Checking capabilities:
// <https://stackoverflow.com/questions/28939652/how-to-detect-sse-sse2-avx-avx2-avx-512-avx-128-fma-kcvi-availability-at-compile>
//
// e.g. "icpc -march=skylake-avx512 -dM -E - < /dev/null | egrep "SSE|AVX" | sort"

// All intel intrinsics:
// <https://software.intel.com/sites/landingpage/IntrinsicsGuide/#>

/*
#ifdef AVX512_4VNNIW
#include <immintrin.h>
#endif 

#include <mmintrin.h>  // MMX technology
#include <xmmintrin.h> // Streaming SIMD intrinsics (SSE)
#include <emmintrin.h> // Streaming SIMD Extensions 2 (SSE2)
*/

// This I took straight from ls1 mardyn's "SIMD_TYPE.h" file, with small name modifications.
#define SIMD_NONE 0
#define SIMD_SSE3 1
#define SIMD_AVX  2
#define SIMD_AVX2 3
#define SIMD_MIC  4

#if defined(__AVX2__) && not defined(__FMA__)//fma should always be existent alongside avx2!!!
  #warn AVX2 enabled, but no FMA found. Please enable fma to use avx2.
#endif

// define symbols for vectorization
#if defined(__MIC__)
  #define SIMD_TYPE SIMD_MIC
#elif defined(__AVX2__) && defined(__FMA__)
  #define SIMD_TYPE SIMD_AVX2
#elif defined(__AVX__) && not defined(AVX128)
  #define SIMD_TYPE SIMD_AVX
#elif defined(__AVX__) && defined(AVX128)
  #define SIMD_TYPE SIMD_SSE3
#elif defined(__SSE3__)
  #define SIMD_TYPE SIMD_SSE3
#else
  #define SIMD_TYPE SIMD_NONE
#endif

#ifdef NOVEC
  #ifdef SIMD_TYPE
    #warn Multiple vectorization methods specified. Will not use vectorization at all!
    #undef SIMD_TYPE
  #endif
  #define SIMD_TYPE SIMD_NONE
#endif

// Include necessary files if we vectorize.
#if SIMD_TYPE==SIMD_AVX or SIMD_TYPE==SIMD_AVX2 or SIMD_TYPE==VCP_VEC_MIC
  #include "immintrin.h"
#elif SIMD_TYPE==SIMD_SSE3
  #include "pmmintrin.h"
#endif

#endif // __SIMD_HPP__GFLOW__