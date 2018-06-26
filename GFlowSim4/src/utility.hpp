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

namespace GFlowSimulation {

  // Define what type of floating point data to use
  typedef float RealType;

  // Unsigned int
  typedef unsigned int uint;

  // Convert an object to a string
  template<typename T> inline string toStr(T obj) {
    stringstream stream;
    stream << obj;
    string str;
    stream >> str;
    return str;
  }

  template<typename T> inline T max(T a, T b) {
    return a>b ? a : b;
  }

  template<typename T> inline T min(T a, T b) {
    return a<b ? a : b;
  }

  template<typename T> inline bool contains(vector<T> vec, T obj) {
    return std::find(vec.begin(), vec.end(), obj) != vec.end();
  }

  /// Get the current time
  inline high_resolution_clock::time_point current_time() {
    return high_resolution_clock::now();
  }

  // high_resolution_clock::time_point
  template<typename T> inline RealType time_span(T end, T start) {
    duration<RealType> span = duration_cast<duration<double> >(end-start);
    return span.count();
  }

  // Print as (hrs):(mins):(sec)
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

  template<typename T> inline void copyArray(const T *source, T *destination, int size) {
    for (int i=0; i<size; ++i) destination[i] = source[i];
  }

} // End namespace GFlowSimulation

// Include this after so RealType is defined
#include "defaultconstants.hpp"
#include "bounds.hpp"

// Include this after so string is defined
#include "exceptions.hpp"

#endif // __UTILITY_HPP__GFLOW__
