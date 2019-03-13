#ifndef __VEC3_HPP__
#define __VEC3_HPP__

namespace GFlowSimulation {

  struct Vec3 {
    Vec3() : x(0), y(0), z(0) {};
    Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {};

    Vec3 operator+(const Vec3 v) const {
      return Vec3(x+v.x, y+v.y, z+v.z);
    }

    Vec3 operator*(float s) const {
      return Vec3(s*x, s*y, s*z);
    }

    float x, y, z;
  };

}
#endif // __VEC3_HPP__