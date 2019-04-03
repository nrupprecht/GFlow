#ifndef __PRINTING_UTILITY_HPP__GFLOW__
#define __PRINTING_UTILITY_HPP__GFLOW__

#include "utility.hpp"
#include <iomanip> // For std::setprecision

namespace GFlowSimulation {

  //! \brief Convert a vector to csv format: x1, x2, x3, ...
  template<typename T> inline string toCSV(const vector<T>& vec, int precision=5) {
    stringstream stream;
    for (int i=0; i<vec.size(); ++i) {
      stream << std::setprecision(precision) << vec[i];
      if (i!=vec.size()-1) stream << std::setprecision(precision) << ",";
    }
    string str;
    stream >> str;
    return str;
  }

  //! \brief Convert a string to an int.
  inline int toInt(string s) {
    stringstream stream;
    int i;
    stream << s;
    stream >> i;
    return i;
  }

  //! \brief Convert a realtype to a string, with set precision.
  inline string toStr(RealType x, int precision) {
    stringstream stream;
    string str;
    stream << std::setprecision(precision) << x;
    stream >> str;
    return str;
  }

  //! \brief Pretty print a string so the decimal point is at a fixed place.
  inline string pprint(RealType x, int lead, int follow) {
    stringstream stream;
    string str;
    stream << x;
    stream >> str;
    // Split string
    bool front = true;
    string leading, following;
    for (auto c : str) {
      if (c=='.') front = false;
      else if (front) leading.push_back(c);
      else following.push_back(c);
    }
    // Work with leading part
    int ds = max(0, static_cast<int>(lead - leading.size()));
    string L(ds, ' ');
    L += leading;
    // Work with following part
    following.resize(follow, ' ');
    // Return the pretty string
    return L+"."+following;
  }

  //! \brief Convert a string to whatever type is specified.
  //!
  //! Use like "double x = convert<double>(str);"
  template<typename T> inline T convert(const string s) {
    stringstream stream;
    T data;
    stream << s;
    stream >> data;
    return data;
  }

  template<typename T> bool writeVectorToFile(vector<T>& vec, string fileName) {
    ofstream fout(fileName);
    if (fout.fail()) return false;
    for (int i=0; i<vec.size(); ++i) {
      fout << vec.at(i);
      if (i!=vec.size()-1) fout << ",";
    }
    return true;
  }

  template<typename T> bool writeVectorVectorToFile(vector<vector<T> >& vec, string fileName) {
    ofstream fout(fileName);
    if (fout.fail()) return false;
    for (int i=0; i<vec.size(); ++i) {
      for (int j=0; j<vec.at(i).size(); ++j) {
        fout << vec.at(i).at(j);
        if (j!=vec.at(i).size()-1) fout << ",";
      }
      fout << endl;
    }
    return true;
  }

  //! \brief Write an array that represents [elements]*[width] data to a .csv file
  bool writeArrayDataToFile(RealType*, int, int, string);

  //! \brief Write a vector to a directory.
  bool writeVectorToDirectory(vector<RealType*>&, const vector<int>&, int, string, const string);

  //! \brief Write a vector to a comma separated string
  string toStrVec(const double*, int);

  //! \brief Write a vector to a comma separated string
  string toStrVec(const float*, int);

  //! \brief Write an integer vector to a comma separated string
  string toStrVec(const int*, int);

  template<typename T> string toStrVec(const T* vec, int number) {
    string str("");
    for (int i=0; i<number; ++i) {
      str += toStr(vec[i]);
      if (i!=number-1) str += ',';
    }
    return str;
  }

}
#endif // __PRINTING_UTILITY_HPP__GFLOW__