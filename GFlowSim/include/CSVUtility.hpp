#ifndef __CSV_UTILITY_HPP__
#define __CSV_UTILITY_HPP__

// Includes
#include "Utility.hpp"
#include "vec2d.hpp"

#include <sys/stat.h> // For linux

namespace GFlow {

  template<typename T> inline string toStr(T d) {
    stringstream stream;
    string str;
    stream << d;
    stream >> str;
    return str;
  }

  inline string toCSV(const int n) { return toStr(n); }

  inline string toCSV(const double n) { return toStr(n); }

  inline string toCSV(const vec2& p) {
    stringstream stream;
    string str;
    stream << p.x << "," << p.y;
    stream >> str;
    return str;
  }

  template<typename... T> string toCSV(const PData& pdata) {
    stringstream stream;
    string str;
    int size = sizeof...(T);
    stream << toCSV(std::get<0>(pdata)) << "," << std::get<1>(pdata) << "," << std::get<2>(pdata) << "," << std::get<3>(pdata) << "," << std::get<4>(pdata) << "," << std::get<5>(pdata);
    stream >> str;
    return str;
  }

  template<typename T> bool printToCSV(string filename, const vector<T>& lst) {
    if (lst.empty()) return true;
    ofstream fout(filename);
    if (fout.fail()) return false;
    for (int i=0; i<lst.size(); i++) fout << toCSV(lst.at(i)) << endl;
    fout.close();
    return true;
  }

  template<typename T> bool printToCSV(string filename, const vector<T>& lst, int iter) {
    stringstream stream;
    stream << filename << iter << ".csv";
    stream >> filename;
    return printToCSV(filename, lst);
  }

  template<typename T> bool printToCSV(string filename, const vector<vector<T> >& lst) {
    if (lst.empty()) return true;
    ofstream fout(filename);
    if (fout.fail()) return false;
    for (int i=0; i<lst.size(); i++) {
      for (int j=0; j<lst.at(i).size(); ++j) {
	fout << lst.at(i).at(j);
	if (j!=lst.at(i).size()-1) fout << ",";
      }
      fout << endl;
    }
    fout.close();
    return true;
  }
  
}

#endif // __CSV_UTILITY_HPP__
