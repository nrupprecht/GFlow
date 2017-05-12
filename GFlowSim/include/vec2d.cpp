#include "vec2d.hpp"

vec2 vec2::operator-(const vec2& v) const {
  return vec2(x-v.x, y-v.y);
}

vec2 vec2::operator+(const vec2& v) const {
  return vec2(x+v.x, y+v.y);
}

RealType vec2::operator*(const vec2& v) const {
  return x*v.x + y*v.y;
}

vec2 operator*(RealType c, const vec2& v) {
  return vec2(c*v.x, c*v.y);
}
