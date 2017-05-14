#ifndef __BOUNDS_HPP__
#define __BOUNDS_HPP__

#include "Utility.hpp"
#include "CSVUtility.hpp"
#include "vec2d.hpp"

struct Bounds {
  Bounds() : left(0), right(0), bottom(0), top(0) {};
  Bounds(RealType l, RealType r, RealType b, RealType t) : left(l), right(r), bottom(b), top(t) {};

  // Data
  RealType left, right, bottom, top;

  // Member and friend functions
  friend std::ostream& operator<<(std::ostream& out, const Bounds& B) {
    out << "{" << B.left << "," << B.right << "," << B.bottom << "," << B.top << "}";
    return out;
  }
  void cut(RealType c) {
    left+=c; right-=c; bottom+=c; top-=c;
  }
  Bounds& operator=(const Bounds& b) {
    left = b.left; right = b.right;
    bottom = b.bottom; top = b.top;
    return *this;
  }
  bool operator==(const Bounds& b) const {
    return b.left==left && b.right==right && b.bottom==bottom && b.top==top;
  }
  bool contains(const vec2 position) const {
    return position.x<right && left<position.x && position.y<top && bottom<position.y;
  }
  bool contains(const vec2 position, const RealType sigma) const {
    return position.x-sigma<right && left<position.x+sigma && position.y-sigma<top && bottom<position.y+sigma;
  }
  bool contains(const RealType x, const RealType y) const {
    return x<right && left<x && y<top && bottom<y;
  }
};

const Bounds NullBounds(0., -1., 0., -1.);

inline string toCSV(const Bounds& b) {
  stringstream stream;
  string str;
  stream << b.left << "," << b.right << "," << b.bottom << "," << b.top;
  stream >> str;
  return str;
}

#endif // __BOUNDS_HPP__
