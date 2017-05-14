#ifndef __VEC_2D_HPP__
#define __VEC_2D_HPP__

#include "Utility.hpp"

using GFlow::RealType;

/*
 * @class vec2
 * Two dimensional vector class
 *
 */
struct vec2 {
  // Default constructor
  vec2() : x(0), y(0) {};
  
  // Initialized constructor
  vec2(RealType x, RealType y) : x(x), y(y) {};
  
  // The actual vector data
  RealType x, y;

  // Operators
  vec2 operator-(const vec2&) const;
  vec2 operator+(const vec2&) const;
  RealType operator*(const vec2&) const;
  friend vec2 operator*(RealType, const vec2&);

  friend std::ostream& operator<<(std::ostream&, const vec2&);
};

const vec2 Zero(0,0);

// Special squaring function for vectors
inline RealType sqr(const vec2& v) { return v*v; }

#endif // __VEC_2D_HPP__
