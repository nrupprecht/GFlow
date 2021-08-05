#ifndef __PRINTING_UTILITY_HPP__GFLOW__
#define __PRINTING_UTILITY_HPP__GFLOW__

#include "utility/utility.hpp"
#include <iomanip> // For std::setprecision

namespace GFlowSimulation {

//! \brief Convert a vector to csv format: x1, x2, x3, ...
template<typename T>
inline std::string toCSV(const std::vector<T> &vec, int precision = 5) {
  std::stringstream stream;
  for (int i = 0; i < vec.size(); ++i) {
    stream << std::setprecision(precision) << vec[i];
    if (i != vec.size() - 1) {
      stream << std::setprecision(precision) << ",";
    }
  }
  string str;
  stream >> str;
  return str;
}

//! \brief Convert a string to an int.
inline int toInt(const std::string &s) {
  std::stringstream stream;
  int i;
  stream << s;
  stream >> i;
  return i;
}

//! \brief Convert a realtype to a string, with set precision.
inline string toStr(RealType x, int precision) {
  std::stringstream stream;
  string str;
  stream << std::setprecision(precision) << x;
  stream >> str;
  return str;
}

//! \brief Pretty print a string so the decimal point is at a fixed place.
inline std::string pprint(RealType x, int lead, int follow) {
  std::stringstream stream;
  std::string str;
  stream << x;
  stream >> str;
  // Split string
  bool front = true;
  std::string leading, following;
  for (auto c : str) {
    if (c == '.') {
      front = false;
    }
    else if (front) {
      leading.push_back(c);
    }
    else {
      following.push_back(c);
    }
  }
  // Work with leading part
  int ds = max(0, static_cast<int>(lead - leading.size()));
  std::string L(ds, ' ');
  L += leading;
  // Work with following part
  following.resize(follow, ' ');
  // Return the pretty string
  return L + "." + following;
}

//! \brief Convert a string to whatever type is specified.
//!
//! Use like "double x = convert<double>(str);"
template<typename T>
inline T convert(const std::string &s) {
  stringstream stream;
  T data(0);
  stream << s;
  stream >> data;
  return data;
}

template<>
inline std::string convert<std::string>(const std::string &s) {
  return s;
}

template<typename T>
bool writeVectorToFile(std::vector<T> &vec, const std::string &fileName) {
  ofstream fout(fileName);
  if (fout.fail()) {
    return false;
  }
  for (int i = 0; i < vec.size(); ++i) {
    fout << vec.at(i);
    if (i != vec.size() - 1) {
      fout << ",";
    }
  }
  return true;
}

template<typename T>
bool writeVectorVectorToFile(std::vector<std::vector<T>> &vec,
                             const std::string& fileName) {
  ofstream fout(fileName);
  if (fout.fail()) {
    return false;
  }
  for (int i = 0; i < vec.size(); ++i) {
    for (int j = 0; j < vec.at(i).size(); ++j) {
      fout << vec.at(i).at(j);
      if (j != vec.at(i).size() - 1) {
        fout << ",";
      }
    }
    fout << endl;
  }
  return true;
}

//! \brief Write an array that represents [elements]*[width] data to a .csv file
bool writeArrayDataToFile(RealType *array,
                          int elements,
                          int width,
                          const std::string& fileName);

//! \brief Write a vector to a directory.
bool writeVectorToDirectory(std::vector<RealType*>& record,
                            const std::vector<int>& elements,
                            int width,
                            std::string dirName,
                            const std::string& fileName);

//! \brief Write a vector to a comma separated string
std::string toStrVec(const double *x, int length);

//! \brief Write a vector to a comma separated string
std::string toStrVec(const float *x, int length);

//! \brief Write an integer vector to a comma separated string
std::string toStrVec(const int *x, int length);

template<typename T>
std::string toStrVec(const T *vec, int number) {
  std::string str;
  for (int i = 0; i < number; ++i) {
    str += toStr(vec[i]);
    if (i != number - 1) {
      str += ',';
    }
  }
  return str;
}

}
#endif // __PRINTING_UTILITY_HPP__GFLOW__
