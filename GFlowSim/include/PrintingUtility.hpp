#ifndef __PRINTING_UTILITY_HPP__
#define __PRINTING_UTILITY_HPP__

#include <vector>
#include <list>

namespace GFlow {

  template<typename T> inline ostream& operator<<(ostream& out, const std::vector<T>& vec) {
    out << "{";
    for (int i=0; i<vec.size(); ++i) {
      out << vec.at(i);
      if (i!=vec.size()-1) out << ",";
    }
    out << "}";
    return out;
  }

  template<typename T> inline ostream& operator<<(ostream& out, const std::list<T>& lst) {
    out << "{";
    int i=0;
    for (const auto v : lst) {
      out << v;
      if (i!=lst.size()-1) out << ",";
    }
    out << "}";
    return out;
  }

  template<typename T, typename U> inline ostream& operator<<(ostream& out, const std::pair<T,U>& vec) {
    out << "{" << vec.first << "," << vec.second << "}";
    return out;
  }

}
#endif // __PRINTING_UTILITY_HPP__
