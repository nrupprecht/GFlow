#ifndef __RAY_TRACE_HPP__
#define __RAY_TRACE_HPP__

#include "../EasyBMP/EasyBMP.h"

#include "vec3.hpp"
#include "kdtree-sphere.hpp"

struct Ray {

  double intersect(Vec3 center, double radius);

  Vec3 origin, orientation;
};

class RayTracer {
public:

  void setBounds(BoundsPack&);

  void addSphere(Vec3, double);

  void addSphere(vector<Vec3>&, vector<double>&);

  const BMP& trace();

private:

  // --- Helper functions
  void partitionSpace();

  BoundsPack bounds;

  //! @brief The image that the tracer creates.
  BMP image;

  //! @brief The position of the camera.
  Vec3 camera;

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