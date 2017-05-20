/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 15, 2017
 *
 */

#ifndef __PRINTING_UTILITY_HPP__
#define __PRINTING_UTILITY_HPP__

// Includes
#include <vector>
#include <list>
#include <ostream>
#include <string>
#include <sstream>

namespace GFlow {

  template<typename S, typename T> inline std::ostream& operator<<(std::ostream& out, const std::pair<S,T> pr) {
    out << "{" << pr.first << "," << pr.second << "}";
    return out;
  }

  template<typename T> inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    out << "{";
    for (int i=0; i<vec.size(); ++i) {
      out << vec.at(i);
      if (i!=vec.size()-1) out << ",";
    }
    out << "}";
    return out;
  }

  template<typename T> inline std::ostream& operator<<(std::ostream& out, const std::list<T>& lst) {
    out << "{";
    int i=0;
    for (const auto v : lst) {
      out << v;
      if (i!=lst.size()-1) out << ",";
    }
    out << "}";
    return out;
  }

  template<typename T, typename U> inline std::ostream& operator<<(std::ostream& out, const std::pair<T,U>& vec) {
    out << "{" << vec.first << "," << vec.second << "}";
    return out;
  }

  // Get rid of all the "e-"s for mathematica, truncate doubles after <precision> decimal places
  inline std::string mmPreproc(std::string s, int precision=-1) {
    int size = s.size();
    std::string out;
    out.reserve(size);
    for (int i=0; i<size; i++) {
      // Get numbers
      if (precision>=0 && isdigit(s.at(i))) {
	std::string number;
	number.push_back(s.at(i));
	i++;
	int count = 0, decimal = 0; // decimal == 1 when we detect a '.'
	for (; i<size && (isdigit(s.at(i)) || s.at(i)=='.'); i++) {
	  if (s.at(i)=='.') { decimal = 1; number.push_back('.'); }
	  else if (decimal==0) number.push_back(s.at(i));
	  else if (count<precision) { number.push_back(s.at(i)); count++; }
	}
	// Done getting the number
	out += number;
	if (i==size) return out; // We are done
	i--;
      }
      // Get nan, replace with 10^10
      else if (s.at(i)=='n' && i+2<size && s.at(i+1)=='a' && s.at(i+2)=='n') {
	out += "10^10";
	i+=2;
      }
      // Get scientific notation
      else if (s.at(i)=='e' && i+2<size && (s.at(i+1)=='-' || s.at(i+1)=='+') && isdigit(s.at(i+2))) out += "*10^";
      else out.push_back(s.at(i));
    }
    return out;
  }
  
  template<typename T> inline std::string mmPreproc(T s, int precision=-1) {
    std::stringstream stream;
    std::string str;
    stream << s;
    stream >> str;
    return mmPreproc(str, precision);
  }

  template<typename T> inline std::string toStr(T x) {
    std::stringstream stream;
    std::string str;
    stream << x;
    stream >> str;
    return str;
  }

  template<typename S, typename T> inline std::string toStr(const std::pair<S,T>& p) {
    std::string first, second;
    std::stringstream stream;
    stream << p.first;
    stream >> first;
    stream.clear();
    stream << p.second;
    stream >> second;
    return ("{"+first+","+second+"}");
  }

}
#endif // __PRINTING_UTILITY_HPP__
