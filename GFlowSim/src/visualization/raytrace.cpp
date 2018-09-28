#include "raytrace.hpp"

double Ray::intersect(Vec3 center, double radius) {

}

void RayTracer::setBounds(BoundsPack& b) {
  bounds = b;
}

void RayTracer::addSphere(Vec3 c, double r) {

}

void RayTracer::addSphere(vector<Vec3>& centers, vector<double>& radii) {

}

void RayTracer::partitionSpace() {

}

const BMP& RayTracer::trace() {



  return image;
}