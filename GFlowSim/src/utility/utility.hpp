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

#include <map>

#include <memory> // For pointers
using std::unique_ptr;
using std::shared_ptr;
using std::weak_ptr;

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

  //! \brief A pair containing a string and an int.
  typedef pair<string, int> SIPair;

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

  template<typename T> inline T max(T a, T b, T c) {
    return max(max(a, b), c);
  }

  template<typename T> inline T max(const vector<T>& v) {
    return *std::max_element(v.begin(), v.end());
  }

  //! Find the min of two objects
  template<typename T> inline T min(T a, T b) {
    return a<b ? a : b;
  }

  template<typename T> inline T min(T a, T b, T c) {
    return min(min(a, b), c);
  }

  template<typename T> inline T min(const vector<T>& v) {
    return *std::min_element(v.begin(), v.end());
  }

  template<typename T> inline T clamp(T a) {
    return a<0 ? T(0) : a;
  }

  template<typename T> inline T sign(T a) {
    return a<0 ? T(-1.) : T(1.);
  }

  //! \brief A helper function that turns a linear address into a tuple address.
  //!
  //! Store in [ [ x00, x01 ... ], [ x10, x11, ...] ... ] -> row major form
  //! This is the template dimension version of the function. It is overloaded for
  //! the zero through three dimensional cases (at the time of this writing).
  inline void getAddress(int linear, int *dims, int *address, int dimensions) {
    // Multiply dims[0] * dims[1] * ... * dims[c-1]
    // If c==0, returns 1
    auto product = [&] (int c) -> int {
      int p=1;
      for (int i=0; i<c; ++i) p *= dims[i];
      return p;
    };

    for (int d=dimensions-1; d>=0; --d) {
      int prod = product(d);
      address[d] = linear / prod;
      linear %= prod;
    }
  }

  inline void getAddressCM(int linear, int *dims, int *address, int dimensions) {
    auto product = [&] (int c) -> int {
      int p = 1;
      for (int i=c; i<dimensions; ++i) p*=dims[i];
      return p;
    };

    for (int d=0; d<dimensions; ++d) {
      int prod = product(d+1);
      address[d] = linear / prod;
      linear %= prod;
    }
  }

  //! Check if a vector contains a specified object
  template<typename T> inline bool contains(const vector<T>& vec, T obj) {
    return std::find(vec.begin(), vec.end(), obj) != vec.end();
  }

  template<typename T> inline bool contains(const std::set<T> s, T obj) {
    return s.find(obj)!=s.end();
  }

  template<typename T, typename S> inline bool contains(const std::map<T, S> m, T key) {
    return m.find(key)!=m.end();
  }

  //! Get the current time
  inline high_resolution_clock::time_point current_time() {
    return high_resolution_clock::now();
  }

  //! Takes two high_resolution_clock::time_point, returns the number of seconds
  template<typename T> inline RealType time_span(T end, T start) {
    duration<double> span = duration_cast<duration<double> >(end-start);
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
    if(seed==0) seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
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
#include "vec.hpp"

// Needs flag from "defaultconstants.hpp"
#if USE_MPI == 1
#include <mpi.h>
#endif

namespace GFlowSimulation {

  //! Computes the volume of a [D]-dimensional sphere - need PI
  inline RealType sphere_volume(const RealType radius, const int D) {
    return pow(PI, D/2.) * pow(radius, D) / tgamma(D/2. + 1.);
  }

  inline RealType inv_sphere_volume(const RealType v, const int D) {
    return pow( tgamma(D/2. + 1.) * v/ pow(PI, D/2.), 1./D);
  }

}

// Include this after so string is defined
#include "exceptions.hpp"

#endif // __UTILITY_HPP__GFLOW__
