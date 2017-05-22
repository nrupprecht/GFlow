/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 19, 2017
 *
 */

#ifndef __READING_UTILITIES_HPP__
#define __READING_UTILITIES_HPP__

#include <vector>
using std::vector;
#include <list>
using std::list;
#include <istream>
using std::istream;
#include "vec2d.hpp"

namespace GFlow {

  // Exception class
  struct BadIstreamRead {
    BadIstreamRead(string t, char b, char g) : type(t), unexpected(b), expected(g) {};
    string type;
    char unexpected, expected;
  };

  inline int getWhiteSpace(istream& in) {
    int count = 1;
    char c;
    in.get(c);
    while(c==' ') {
      in.get(c);
      ++count;
    }
    in.putback(c);
    --count;
    return count;
  }

  inline istream& operator>>(istream& in, vec2& vec) {
    // Storage
    char c;
    RealType X, Y;
    // Get the vector data
    in.get(c);
    getWhiteSpace(in);
    if (c!='{') throw BadIstreamRead("vec2", c, '{');
    getWhiteSpace(in);
    in >> X;
    in >> c; 
    if (c!=',') throw BadIstreamRead("vec2", c, ',');
    in >> Y;
    getWhiteSpace(in);
    in.get(c);
    if (c!='}') throw BadIstreamRead("vec2", c, '}');
    // Set vec
    vec.x = X; vec.y = Y;
    // Return the istream
    return in;
  }

  template<typename T> inline istream& operator>>(istream& in, vector<T>& vec) {
    // Clear the vector
    vec.clear();
    // Storage
    char c;
    T data;
    bool reading = true;
    // Get the vector data    
    in.get(c);
    while (!in.eof() && reading) {
      if (c==' ');
      else if (c==',');
      else if (isdigit(c)) {
	in.putback(c);
	in >> data;
	vec.push_back(data);
      }
      else if (c=='}') {
	reading = false;
	continue;
      }
      else throw BadIstreamRead("vec2", c, '*')
      in.get(c);
    }
    // Return the istream
    return in;
  }  
  
  inline string printChar(char c) {
    if (c=='\n') return "\\n";
    if (c=='\r') return "\\r";
    if (c=='\t') return "\\t";
    string str;
    str.push_back(c);
    return str;
  }

}
#endif // __READING_UTILITIES_HPP__
