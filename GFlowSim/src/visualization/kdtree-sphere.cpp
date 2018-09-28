#include "kdtree-sphere.hpp"

KDTreeNode::~KDTreeNode() {
  if (min) delete min;
  if (max) delete max;
}

KDTreeSphere::~KDTreeSphere() {

}

void KDTreeSphere::construct(const vector<Vec3>& centers, const vector<double>& radii) {

}

inline void KDTreeSphere::constructByHalving() {
  
}