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

#include <ctime>

#include <cmath>

namespace GFlowSimulation {

  // Define what type of floating point data to use
  typedef float RealType;

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

}

// Include this after so RealType is defined
#include "defaultconstants.hpp"
#include "bounds.hpp"
#include "printingutility.hpp"

// Include this after so string is defined
#include "exceptions.hpp"

#endif // __UTILITY_HPP__GFLOW__