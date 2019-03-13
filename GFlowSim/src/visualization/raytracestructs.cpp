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
  }

  Sphere::Sphere(const float *c, const float r) : radius(r) {
    copyVec(c, center, 3);
  }

  Sphere::Sphere(const Sphere& s) : radius(s.radius) {
    copyVec(s.center, center, 3);
  }

  Sphere& Sphere::operator=(const Sphere& s) {
    radius = s.radius;
    copyVec(s.center, center, 3);
    return *this;
  }

  bool Sphere::intersect(const Ray& ray, float* point, float& distance_near, float& distance_far) const {
    // Vector from ray origin to sphere center: E = sphere.center - ray.origin
    float E[3];
    subtractVec(center, ray.origin, E, 3);
    // Distance from ray origin to point of closest approach.
    float v = dotVec(ray.orientation, E, 3);
    // If v<0, the point of closest approach is behind you, and there is no intersection.
    if (v<=0) return false;
    // Find the distance from the point of closest approach to the intersection point.
    float d = sqr(radius) - (dotVec(E, E, 3)-sqr(v));
    if (d<0) return false;
    d = sqrtf(d);

    // Just find the closer point. This assumes the ray does not start within the sphere.
    distance_near = v - d;
    distance_far  = v + d; 
    scalarMultVec(distance_near, ray.orientation, point, 3);
    plusEqVec(point, ray.origin, 3);
    // There was an intersection.
    return true;
  }

  void swap(Sphere& s1, Sphere& s2) {
    std::swap(s1.center[0], s2.center[0]);
    std::swap(s1.center[1], s2.center[1]);
    std::swap(s1.center[2], s2.center[2]);
    std::swap(s1.radius, s2.radius);
  }

}