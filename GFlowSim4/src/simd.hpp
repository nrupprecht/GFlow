#pragma once

#ifndef __SIMD_HPP__GFLOW__
#define __SIMD_HPP__GFLOW__


// Checking capabilities:
// <https://stackoverflow.com/questions/28939652/how-to-detect-sse-sse2-avx-avx2-avx-512-avx-128-fma-kcvi-availability-at-compile>
//
// e.g. "icpc -march=skylake-avx512 -dM -E - < /dev/null | egrep "SSE|AVX" | sort"

// All intel intrinsics:
// <https://software.intel.com/sites/landingpage/IntrinsicsGuide/#>

#ifdef AVX512_4VNNIW
#include <immintrin.h>
#endif 

#include <mmintrin.h>  // MMX technology
#include <xmmintrin.h> // Streaming SIMD intrinsics (SSE)
#include <emmintrin.h> // Streaming SIMD Extensions 2 (SSE2)

#if SIMD_TYPE==NO_SIMD

#endif

#endif // __SIMD_HPP__GFLOW__