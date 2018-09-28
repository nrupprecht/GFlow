#ifndef __VEC3_HPP__
#define __VEC3_HPP__

struct CameraMatrix {

  
  
  double *array;
};

struct Vec3 {
  Vec3() : x(0), y(0), z(0) {};
  Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};

  Vec3 operator+(const Vec3 v) const {
    return Vec3(x+v.x, y+v.y, z+v.z);
  }

  Vec3 operator*(double s) const {
    return Vec3(s*x, s*y, s*z);
  }

  double x, y, z;
};

#endif // __VEC3_HPP__