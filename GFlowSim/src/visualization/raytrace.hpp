#ifndef __RAY_TRACE_HPP__
#define __RAY_TRACE_HPP__

#include "../EasyBMP/EasyBMP.h"

#include "vec3.hpp"
#include "kdtree-sphere.hpp"

class RayTracer {
public:

  void addSphere(Vec3, double);

  const BMP& trace();

private:

  BMP image;

  //! @brief The position of the camera, or "eye."
  Vec3 eye;

  //! @brief The direction of maximum light
  Vec3 lightDirection;

  //! @brief surfaces with normals within this angle of the light direction "directly reflect" the light source.
  double reflection_angle;

  //! @brief A vector of the centers of all the spheres.
  vector<Vec3> centers;
  //! @brief A vector of the radii of all the spheres.
  vector<double> radii;
};

#endif