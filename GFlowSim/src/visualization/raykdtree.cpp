#include "raykdtree.hpp"
// Other files
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  RayKDTreeNode::~RayKDTreeNode() {
    if (lower)  delete lower;
    lower = nullptr;
    if (higher) delete higher;
    higher = nullptr;
  }

  void RayKDTreeNode::insert(Sphere *sphere) {
    // If this is a leaf node, insert the sphere
    if (lower==nullptr) {
      sphere_list.push_back(sphere);
      return;
    }
    // Else, check which children the sphere should be inserted into.
    if (sphere->center[split_dim] - sphere->radius < split_value) lower->insert(sphere);
    if (split_value < sphere->center[split_dim] + sphere->radius) higher->insert(sphere);
  }

  void RayKDTreeNode::empty() {
    sphere_list.clear();
    if (lower) lower->empty();
    if (higher) higher->empty();
  }

  Sphere* RayKDTreeNode::traverse(const Ray& ray, float *point, float& distance_near, float& distance_far, float tmin, float tmax) const {
    // If this is a terminal node, search through all spheres for an intersection.
    if (lower==nullptr) {
      float distance = 10000, d_far = 10000;
      distance_near = distance_far = 10000;
      float test_point[3];
      Sphere *min_sphere = nullptr;
      // Go through all spheres in this node.
      for (const auto& sphere : sphere_list) {
        if (sphere->intersect(ray, test_point, distance, d_far, tmin) && distance<distance_near) {
          min_sphere = sphere;
          distance_near = distance;
          distance_far  = d_far;
          copyVec(test_point, point, 3);
        }
      }
      // If we found nothing, this just returns nullptr.
      return min_sphere;
    }

    // Else, march the ray through any child nodes that it intersects with.
    float v = ray.orientation[split_dim];
    if (v!=0) { // Make sure ray is not parallel to splitting plane
      // Calculate the t value for which the ray intersects with the splitting plane.
      float t_split = (split_value - ray.origin[split_dim])/v;

      // In eiher case, we will need to know which child domain the ray enters first.
      float X[3];
      scalarMultVec(tmin, ray.orientation, X, 3);
      plusEqVec(X, ray.origin, 3);
      bool lower_first = (X[split_dim] < split_value);

      // If t_split > tmax, the ray exits the bounding box before entering the other child domain.
      // If t_split < tmin, the ray is going the wrong direction to go into one of the child domains.
      // In either case, we only need to check with one child domain. We can test which one is relevant by seeing which child 
      // domain the entrance point is in.
      if (t_split>tmax || t_split<tmin) {
        // Check which side of the splitting coordinate X is on
        if (lower_first) return lower->traverse(ray, point, distance_near, distance_far, tmin, tmax);
        else return higher->traverse(ray, point, distance_near, distance_far, tmin, tmax);
      }
      
      // Otherwise, the ray intersects with both child domains. The only thing left to do is to figure out which domain comes first.
      // We can do this by seeing which child domain the entrance point is in
      Sphere *sphere = nullptr;
      if (lower_first) {
        sphere = lower->traverse(ray, point, distance_near, distance_far, tmin, tmax);
        // If nothing was found in the first child, search the second child.
        if (sphere==nullptr) return higher->traverse(ray, point, distance_near, distance_far, tmin, tmax);
        // If something was found return it.
        else return sphere;
      }
      else {
        sphere = higher->traverse(ray, point, distance_near, distance_far, tmin, tmax);
        // If nothing was found in the first child, search the second child.
        if (sphere==nullptr) return lower->traverse(ray, point, distance_near, distance_far, tmin, tmax);
        // If something was found return it.
        else return sphere;
      }

    }
    else {
      // Only intersects the lower node
      if (ray.origin[split_dim] < split_value) return lower->traverse(ray, point, distance_near, distance_far, tmin, tmax);
      // Only intersects the higher node
      else return higher->traverse(ray, point, distance_near, distance_far, tmin, tmax);
    }
  }

  int RayKDTreeNode::getMaxDepth() const {
    if (lower==nullptr) return 1;
    return 1 + max(lower->getMaxDepth(), higher->getMaxDepth());
  }

  int RayKDTreeNode::getTotalSpheres() const {
    if (lower==nullptr) return sphere_list.size();
    return lower->getTotalSpheres() + higher->getTotalSpheres();
  }

  RayKDTree::~RayKDTree() {
    clear();
  }

  void RayKDTree::insert(Sphere *sphere) {
    if (head) head->insert(sphere);
  }

  Sphere* RayKDTree::traverse(const Ray& ray, float *point, float& distance_near, float& distance_far, bool& intersect, float& tmin, float& tmax) const {
    // If the tree is empty, return nullptr.
    if (head==nullptr) return nullptr;
    // Find the tmin and tmax of the ray through the system bounding box.
    tmin = 0; tmax = 0;
    // If the ray does not intersect with the bounding box, set the flag and return.
    if (!getRayIntersectionParameters(ray, tmin, tmax)) {
      intersect = false;
      return nullptr;
    }
    // The ray did intersect with the bounding box.
    intersect = true;
    // Traverse the data structure.
    return head->traverse(ray, point, distance_near, distance_far, tmin, tmax);
  }

  bool RayKDTree::getRayIntersectionParameters(const Ray& ray, float& tmin, float& tmax) const {
    float tmin_d[3], tmax_d[3];
    // Find tmin, tmax for each dimension.
    for (int i=0; i<3; ++i) {
      if (ray.orientation[i]==0) {
        // No intersection.
        if (ray.origin[i]<=bounds.min[i] || bounds.max[i]<=ray.origin[i]) return false;
        tmin_d[i] = 0;
        tmax_d[i] = 10000; // Use as "infinity."
      }
      else if (ray.orientation[i]>0) {
        tmin_d[i] = max(static_cast<float>((bounds.min[i] - ray.origin[i])/ray.orientation[i]), 0.f);
        tmax_d[i] = max(static_cast<float>((bounds.max[i] - ray.origin[i])/ray.orientation[i]), 0.f);
      }
      else { // ray.orientation[i]<0
        tmin_d[i] = max(static_cast<float>((bounds.max[i] - ray.origin[i])/ray.orientation[i]), 0.f);
        tmax_d[i] = max(static_cast<float>((bounds.min[i] - ray.origin[i])/ray.orientation[i]), 0.f);
      }
    }
    // Find actual tmin, tmax
    tmin = max(tmin_d[0], tmin_d[1], tmin_d[2]);
    tmax = min(tmax_d[0], tmax_d[1], tmax_d[2]);
    // If the ray does not intersect with the bounding box, set the flag and return.
    if (tmax<=tmin) return false;
    // The ray did intersect with the bounding box.
    return true;
  }

  void RayKDTree::clear() {
    if (head) delete head;
    head = nullptr;
  }

  void RayKDTree::empty() {
    if (head) head->empty();
  }

  bool RayKDTree::isCreated() {
    return head!=nullptr;
  }

  int RayKDTree::getMaxDepth() const {
    if (head) return head->getMaxDepth();
    else return 0;
  }

  int RayKDTree::getTotalSpheres() const {
    if (head) return head->getTotalSpheres();
    else return 0;
  }

}