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

  // PData records all the data you need to print particles:
  // { pos x, pos y, sigma, theta, interaction, 'color' }
  typedef std::tuple<double, double, double, double, double, double> PData;
}

#endif // __UTILITY_HPP__
