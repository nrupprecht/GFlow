#ifndef __PRINTING_UTILITY_HPP__GFLOW__
#define __PRINTING_UTILITY_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  template<typename T, int dimensions>
  string toStr(const vec<T, dimensions> v) {
    stringstream stream;
    stream << "{";
    for (int d=0; d<dimensions; ++d) {
      stream << v[d];
      if (d!=dimensions-1) stream << ",";
    }
    stream << "}";
    string str;
    stream >> str;
    return str;
  }

  template<typename T> inline string toCSV(const vector<T>& vec) {
    stringstream stream;
    for (int i=0; i<vec.size(); ++i) {
      stream << vec[i];
      if (i!=vec.size()-1) stream << ",";
    }
    string str;
    stream >> str;
    return str;
  }

  //! @todo Turn all these functions into free-standing (inline) functions.
  struct PrintingUtility {

    template<typename T> static bool writeVectorToFile(vector<T>& vec, string fileName) {
      ofstream fout(fileName);
      if (fout.fail()) return false;
      for (int i=0; i<vec.size(); ++i) {
        fout << vec.at(i);
        if (i!=vec.size()-1) fout << ",";
      }
      return true;
    }

    template<typename T> static bool writeVectorVectorToFile(vector<vector<T> >& vec, string fileName) {
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

    // Write an array that represents [elements]*[width] data to a .csv file
    static bool writeArrayDataToFile(RealType*, int, int, string);

    // [dirName] is the name of the directory that we should create our new directory of data in
    // [fileName] is the name of the new directory, and the files in that directory are called [fileName][#].csv
    static bool writeVectorToDirectory(vector<RealType*>&, const vector<int>&, int, string, const string);

    static bool writeVerletListToDirectory(const class VerletList&, const string);

    // Write a vector to a comma separated string
    static string toStrVec(const RealType*);

    // Write an integer vector to a comma separated string
    static string toStrVec(const int*);

    template<typename T> static string toStrVec(const T* vec, int number) {
      string str("");
      for (int i=0; i<number; ++i) {
        str += toStr(vec[i]);
        if (i!=number-1) str += ',';
      }
      return str;
    }

  };

}
#endif // __PRINTING_UTILITY_HPP__GFLOW__