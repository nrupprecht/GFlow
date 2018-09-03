#ifndef __UTILITY_HPP__GFLOW__
#define __UTILITY_HPP__GFLOW__

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <vector>
using std::vector;
using std::pair;

#include <list>
using std::list;

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

#include <sys/stat.h> // For mkdir

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include <ctime>

#include <cmath>

#include <algorithm>

#include <random>

#include <stdlib.h> // For aligned_alloc

#include <algorithm> // For std::copy, std::fill, etc.

#include <set>

#include "simd_types.hpp"

/**
*  \namespace GFlowSimulation The GFlow simulation namespace
*
*  This namespace encapsulated all the objects, typedefs, etc. used by the
*  GFlow simulation. Almost every object is contained within this namespace.
*/
namespace GFlowSimulation {

  // Define what type of floating point data to use
  typedef float RealType;

  // Unsigned int
  typedef unsigned int uint;

  // --- Pair typedefs

  //! Real pair
  typedef pair<RealType, RealType> RPair;

  //! Int pair
  typedef pair<int, int> IPair;

  //! Real-Int pair
  typedef pair<RealType, int> RIPair;

  //! Int-Real pair
  typedef pair<int, RealType> IRPair;

  // --- Common functions

  //! Convert an object to a string
  template<typename T> inline string toStr(T obj) {
    stringstream stream;
    stream << obj;
    string str;
    stream >> str;
    return str;
  }

  //! Find the max of two objects
  template<typename T> inline T max(T a, T b) {
    return a>b ? a : b;
  }

  //! Find the min of two objects
  template<typename T> inline T min(T a, T b) {
    return a<b ? a : b;
  }

  //! Check if a vector contains a specified object
  template<typename T> inline bool contains(const vector<T>& vec, T obj) {
    return std::find(vec.begin(), vec.end(), obj) != vec.end();
  }

  template<typename T> inline bool contains(const std::set<T> s, T obj) {
    return s.find(obj)!=s.end();
  }

  //! Get the current time
  inline high_resolution_clock::time_point current_time() {
    return high_resolution_clock::now();
  }

  //! Takes two high_resolution_clock::time_point, returns the number of seconds
  template<typename T> inline RealType time_span(T end, T start) {
    duration<RealType> span = duration_cast<duration<double> >(end-start);
    return span.count();
  }

  //! Print a time given in seconds as (hrs):(mins):(sec)
  inline string printAsTime(double seconds) {
    stringstream stream;
    string str;
    // Just print seconds as usual
    if (seconds<60) {
      stream << seconds;
      stream >> str;
      return str;
    }
    // Print full
    int hours = floor(seconds/3600.);
    seconds -= 3600*hours;
    int minutes = floor(seconds/60.);
    seconds -= 60*minutes;
    // Round seconds
    int sec = seconds;
    // Print hours
    if (hours==0);
    else if (hours<10) stream << "0" << hours << ":";
    else stream << hours << ":";
    // Print minutes
    if (minutes<10 && hours>0) stream << "0" << minutes << ":";
    else if (minutes<10) stream << minutes << ":";
    else stream << minutes << ":";
    // Print seconds
    if (seconds<10) stream << "0" << sec;
    else stream << sec;
    stream >> str;
    return str;
  }

  //! Copy an array of size [size] from [source] to [destination]
  template<typename T> inline void copyArray(const T *source, T *destination, int size) {
    //  for (int i=0; i<size; ++i) destination[i] = source[i];
    std::copy(source, source+size, destination);
  }

  //! Global random number generator
  static std::mt19937 global_generator;
  //! Global normal distribution
  static std::normal_distribution<RealType> global_normal_dist(0., 1.);

  inline void seedNormalDistribution(unsigned seed=0) {
    if(seed==0) seed = std::chrono::system_clock::now().time_since_epoch().count();
    global_generator = std::mt19937(seed);
  }

  //! Returns a random normal double from the global normal distribution
  inline double randNormal() {
    return global_normal_dist(global_generator);
  }

} // End namespace GFlowSimulation

// Include this after so RealType is defined
#include "defaultconstants.hpp"
#include "bounds.hpp"

// Needs flag from "defaultconstants.hpp"
#if USE_MPI == 1
#include <mpi.h>
#endif

namespace GFlowSimulation {

  //! Computes the volume of a [D]-dimensional sphere - need PI
  inline RealType sphere_volume(const RealType radius, const int D=DIMENSIONS) {
    return pow(PI, D/2.) * pow(radius, D) / tgamma(D/2. + 1.);
  }

}

// Include this after so string is defined
#include "exceptions.hpp"

// Other
#include "vector_array.hpp"

#endif // __UTILITY_HPP__GFLOW__
