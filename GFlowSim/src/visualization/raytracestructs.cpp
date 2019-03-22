#include "raytracestructs.hpp"
// Other files
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  Ray::Ray(const float *o, const float *n) {
    copyVec(o, origin, 3);
    copyVec(n, orientation, 3);
  }

  Sphere::Sphere() : radius(1.f) {
    zeroVec(center, 3);
    // Set color reflectivity
    color_reflectivity[0] = color_reflectivity[1] = color_reflectivity[2] = 1.; // White spheres by default
  }

  Sphere::Sphere(const float *c, const float r) : radius(r) {
    copyVec(c, center, 3);
    // Set color reflectivity
    color_reflectivity[0] = color_reflectivity[1] = color_reflectivity[2] = 1.; // White spheres by default
  }

  Sphere::Sphere(const Sphere& s) : radius(s.radius) {
    copyVec(s.center, center, 3);
    // Set color reflectivity
    copyVec(s.color_reflectivity, color_reflectivity, 3);
  }

  Sphere& Sphere::operator=(const Sphere& s) {
    radius = s.radius;
    // Set sphere center
    copyVec(s.center, center, 3);
    // Set color reflectivity
    copyVec(s.color_reflectivity, color_reflectivity, 3);
    // Return the sphere
    return *this;
  }

  bool Sphere::intersect(const Ray& ray, float* point, float& distance_near, float& distance_far, const float tmin) const {
    // Vector from ray origin to sphere center: E = sphere.center - ray.origin
    float E[3];
    E[0] = center[0] - ray.origin[0];
    E[1] = center[1] - ray.origin[1];
    E[2] = center[2] - ray.origin[2];
    // Distance from ray origin to point of closest approach.
    float v = ray.orientation[0]*E[0] + ray.orientation[1]*E[1] + ray.orientation[2]*E[2];
    // If v<0, the point of closest approach is behind you, and there is no intersection.
    if (v<=0) return false;
    // Find the distance from the point of closest approach to the intersection point.
    float d = sqr(radius) - (E[0]*E[0] + E[1]*E[1] + E[2]*E[2] - v*v);
    if (d<0) return false;
    d = sqrtf(d);

    // Just find the closer point. This assumes the ray does not start within the sphere.
    distance_near = v - d;
    distance_far  = v + d; 

    // Check if this is a real intersection, i.e. part of ray intersects with the sphere while inside the scene bounding box.
    if (distance_far<tmin) return false;

    // Find the point of intersection
    point[0] = distance_near*ray.orientation[0] + ray.origin[0];
    point[1] = distance_near*ray.orientation[1] + ray.origin[1];
    point[2] = distance_near*ray.orientation[2] + ray.origin[2];

    // There was an intersection.
    return true;
  }

  void swap(Sphere& s1, Sphere& s2) {
    // Swap center and color reflectivity
    for (int i=0; i<3; ++i) {
      std::swap(s1.center[i], s2.center[i]);
      std::swap(s1.color_reflectivity[i], s2.color_reflectivity[i]);
    }
    // Swap radii
    std::swap(s1.radius, s2.radius);
  }

}