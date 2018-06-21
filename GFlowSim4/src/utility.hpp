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

namespace GFlowSimulation {

  // Define what type of floating point data to use
  typedef float RealType;

}

// Include this after so RealType is defined
#include "defaultconstants.hpp"
#include "bounds.hpp"

#endif // __UTILITY_HPP__GFLOW__