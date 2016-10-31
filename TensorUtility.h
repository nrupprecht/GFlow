#ifndef TENSOR_UTILITY_H
#define TENSOR_UTILITY_H

#include "Utility.h"

/// For use with the CRC
#include"/afs/crc.nd.edu/x86_64_linux/intel/15.0/mkl/include/mkl.h"

#include <math.h>   // For sqrt
#include <stdlib.h> // For drand48
#include <ctime>    // For timing

#include <vector>   // For vector
using std::vector;
using std::pair;

#include <string>   // For string
using std::string;

#include <sstream>  // For string stream
using std::stringstream;

#include <ostream>
using std::ostream;
#include <fstream>  // For file i/o

// For debugging
#include <iostream>
using std::cout;
using std::endl;

// Common function template
typedef double (*function) (double);

template<typename T> string print(const vector<T>& array) {
  if (array.size()==0) return "{}";
  string str;
  stringstream stream;
  stream << "{";
  for (int i=0; i<array.size()-1; i++) stream << array.at(i) << ",";
  stream << array.at(array.size()-1) << "}";
  stream >> str;
  return str;
}

template<typename T, typename U> ostream& operator<<(ostream& out, pair<T, U> P) {
  out << "{" << P.first << "," << P.second << "}";
  return out;
}

#endif
