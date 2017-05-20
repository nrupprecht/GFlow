/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 12, 2017
 *
 */

#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__

// Use MPI
#define USE_MPI 1

// Includes  
#include <iostream>
using std::cout;
using std::endl;

#include <ostream>
using std::ostream;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::stringstream;

#include <vector>
using std::vector;
using std::pair;

#include <list>
using std::list;
  
#include <tuple>
using std::tuple;

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include <string>
using std::string;

#include <random>

#include <math.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace GFlow {

  // What type to use as our real type
  typedef double RealType;

  // What to use for Verlet list
  typedef list<int>          VListSubType;
  typedef list<VListSubType> VListType;

  // What to use for Wall list
  typedef list<int>          WListSubType;
  typedef list<WListSubType> WListType;

  // Constants
  const double PI = 3.14159265358979;

  // Squaring function
  template<typename T> inline T sqr(const T& value) {
    return value*value;
  }

  // Min function
  template<typename T> inline T min(const T& a, const T& b) {
    return a<b ? a : b;
  }

  // Max function
  template<typename T> inline T max(const T& a, const T& b) {
    return a<b ? b : a;
  }

  // Clamp function
  inline RealType clamp(const RealType x) {
    return x<0 ? 0 : x;
  }

  inline RealType sign(const RealType x) {
    return x<0 ? -1. : 1.;
  }

  /// Get the current time
  inline auto current_time() {
    return high_resolution_clock::now();
  }

  // high_resolution_clock::time_point
  template<typename T> inline RealType time_span(T end, T start) {
    duration<RealType> span = duration_cast<duration<double>>(end-start);
    return span.count();
  }

  // PData records all the data you need to print particles:
  // { pos x, pos y, sigma, theta, interaction, 'color' }
  typedef std::tuple<RealType, RealType, RealType, RealType, RealType, RealType> PData;

  // Random number generators
  static std::mt19937 generator;
  static std::normal_distribution<RealType> normal_dist(0., 1.);
  static std::poisson_distribution<int> poisson_dist();
}

#endif // __UTILITY_HPP__
