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

  // Get rid of all the "e-"s for mathematica, truncate doubles after <precision> decimal places
  inline string mmPreproc(string s, int precision=-1) {
    int size = s.size();
    string out;
    out.reserve(size);
    for (int i=0; i<size; i++) {
      // Get numbers
      if (precision>=0 && isdigit(s.at(i))) {
	string number;
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
  
  template<typename T> inline string mmPreproc(T s, int precision=-1) {
    stringstream stream;
    string str;
    stream << s;
    stream >> str;
    return mmPreproc(str, precision);
  }

  template<typename T> inline string toStr(T x) {
    stringstream stream;
    string str;
    stream << x;
    stream >> str;
    return str;
  }

}
#endif // __PRINTING_UTILITY_HPP__
