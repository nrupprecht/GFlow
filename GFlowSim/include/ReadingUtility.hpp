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
#include <map>
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

  inline void getWhiteSpace(istream& in) {
    char c;
    in.get(c);
    while(!in.eof() && (c==' ' || c=='\n' || c=='\r')) in.get(c);
    in.putback(c);
  }

  inline void getNext(istream& in, char& c) {
    getWhiteSpace(in);
    in.get(c);
  }

  inline istream& operator>>(istream& in, vec2& vec) {
    // Storage
    char c;
    RealType X, Y;
    // Get the vector data
    getNext(in, c);
    if (c!='{') throw BadIstreamRead("vec2", c, '{');
    getWhiteSpace(in);
    in >> X;
    getNext(in,c);
    if (c!=',') throw BadIstreamRead("vec2", c, ',');
    getWhiteSpace(in);
    in >> Y;
    getNext(in,c);
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
    getNext(in,c);
    // Get opening '{'
    if (c!='{') throw BadIstreamRead("vec2", c, '{');
    getNext(in,c);
    while (!in.eof() && reading) {
      if (c==' ');
      else if (c==',');
      else if (c=='}') {
	reading = false;
	continue;
      }
      else {
	in.putback(c);
	in >> data;
	vec.push_back(data);
      }
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

  template<typename T> void getValue(istream& in, T& data, const std::map<string, string> variables) {
    string str;
    in >> str;
    auto var = variables.find(str);
    if (var!=variables.end()) { // This is the keyword for a variable
      stringstream stream;
      stream << var->second;
      stream >> data;
    }
    else {
      stringstream stream;
      stream << str;
      stream >> data;
    }
  }

}
#endif // __READING_UTILITIES_HPP__
