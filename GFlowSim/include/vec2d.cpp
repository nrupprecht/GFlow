#include "vec2d.hpp"

vec2 vec2::operator-(const vec2& v) const {
  return vec2(x-v.x, y-v.y);
}

vec2& vec2::operator-=(const vec2& v) {
  x -= v.x; y -= v.y;
  return *this;
}

vec2 vec2::operator+(const vec2& v) const {
  return vec2(x+v.x, y+v.y);
}

vec2& vec2::operator+=(const vec2& v) {
  x += v.x; y += v.y;
  return *this;
}

RealType vec2::operator*(const vec2& v) const {
  return x*v.x + y*v.y;
}

vec2 operator*(const RealType c, const vec2& v) {
  return vec2(c*v.x, c*v.y);
}

std::ostream& operator<<(std::ostream& out, const vec2& v) {
  out << "{" << v.x << "," << v.y << "}";
  return out;
}
