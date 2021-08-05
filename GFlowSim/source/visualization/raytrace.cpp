#include "visualization/raytrace.hpp"
// Other files.
#include "utility/vectormath.hpp"
#include "visualization/palette.hpp" // For color definitions

using namespace GFlowSimulation;

RayTrace::RayTrace()
    : pix_x(768),
      pix_y(768),
      target_leaf_size(50),
      min_dimension(0.f),
      use_tree(true),
      adaptive_split(true),
      chop_spheres(false) {
  // Camera center
  camera_center[0] = camera_center[1] = camera_center[2] = 0;
  // Camera orientation
  camera_orientation[0] = 1;
  camera_orientation[1] = camera_orientation[2] = 0;
  // Camera right
  camera_up[0] = camera_up[2] = 0;
  camera_up[1] = 1;
  // Create light direction
  light_direction[0] = 1. / sqrt(3);
  light_direction[1] = 1. / sqrt(3);
  light_direction[2] = 1. / sqrt(3);
  normalizeVec(light_direction, 3);
  // Light color - white by default
  light_color[0] = light_color[1] = light_color[2] = 1.f;
}

void RayTrace::initialize() {
  // Make sure camera vectors are good.
  normalizeVec(camera_orientation, 3);
  normalizeVec(camera_up, 3);
  // Correct camera right vector
  float v = dotVec(camera_up, camera_orientation, 3);
  if (v != 0) {
    if (v == 1) { // camera_up == camera_orientation
      camera_up[0] += 1; // Now they will be different. (Hopefully)
      normalizeVec(camera_up, 3);
      v = dotVec(camera_up, camera_orientation, 3);
    }
    float proj[3];
    scalarMultVec(v, camera_orientation, proj, 3);
    minusEqVec(camera_up, proj, 3);
    normalizeVec(camera_up, 3);
  }

  // Set up KD tree. If a tree already exists, just add the spheres.
  createKDTree();
}

void RayTrace::render() {
  // Calculate camera up
  float camera_right[3];
  crossVec(camera_orientation, camera_up, camera_right);
  // Set bitmap size
  image.SetSize(pix_x, pix_y);
  // Timing
  auto start_time = current_time();
  // Pixel for storing color
  RGBApixel color;
  // Ray orientation, ray direction vectors, accumulator
  float orientation[3], vx[3], vy[3], acc[3];
  // Ray trace each pixel
  for (int y = 0; y < pix_y; ++y) {
    for (int x = 0; x < pix_x; ++x) {
      // Calculate pixel shift
      float dx = (float(x) - pix_x / 2.f) / float(pix_x);
      float dy = (float(y) - pix_y / 2.f) / float(pix_y);
      // Calculate pixel vector
      scalarMultVec(dx, camera_right, vx, 3);
      scalarMultVec(dy, camera_up, vy, 3);
      // Calculate screen position
      addVec(camera_orientation, vx, acc, 3);
      addVec(acc, vy, orientation, 3);
      // Normalize the camera vector
      normalizeVec(orientation, 3);
      // Create the ray
      Ray ray(camera_center, orientation);
      // Trace the ray to get the color
      color = trace(ray);
      // Set the color - because of how image is laid out, Y = pix_y - y - 1
      image.SetPixel(x, pix_y - y - 1, color);
    }
  }
  // Timing
  auto end_time = current_time();
  render_time += time_span(end_time, start_time);
}

void RayTrace::saveImage(std::string fileName) {
  image.WriteToFile(fileName.c_str());
}

void RayTrace::reserve(int capacity) {
  sphere_list.reserve(capacity);
}

void RayTrace::addSphere(const float *center, const float sigma) {
  sphere_list.push_back(Sphere(center, sigma));
}

void RayTrace::addSphere(const float *center, const float sigma, const RGBApixel color) {
  Sphere sphere(center, sigma);
  sphere.color_reflectivity[0] = color.Red / 255.;
  sphere.color_reflectivity[1] = color.Green / 255.;
  sphere.color_reflectivity[2] = color.Blue / 255.;
  sphere_list.push_back(sphere);
}

void RayTrace::setCameraPlacement(const float *c) {
  copyVec(c, camera_center, 3);
}

void RayTrace::setCameraOrientation(const float *o) {
  copyVec(o, camera_orientation, 3);
}

void RayTrace::setBounds(const Bounds &bnds) {
  kdTree.bounds = bnds;
}

void RayTrace::setBoundsScaled(const Bounds &bnds, real scale) {
  kdTree.bounds = bnds;
  real width[3] = {bnds.wd(0), bnds.wd(1), bnds.wd(2)};
  for (int d = 0; d < 3; ++d) {
    kdTree.bounds.min[d] -= 0.5 * scale * width[d];
    kdTree.bounds.max[d] += 0.5 * scale * width[d];
  }
}

void RayTrace::setResolution(int res) {
  pix_x = pix_y = res;
}

void RayTrace::empty() {
  kdTree.empty();
  sphere_list.clear();
}

double RayTrace::getRenderTime() {
  return render_time;
}

RGBApixel RayTrace::trace(const Ray &ray) const {
  // Distance to intersection point
  float distance_near = 0, distance_far = 0;
  // Intersection point
  float point[3];
  // Did the ray intersect with the scene bounding box?
  bool intersect = false;
  // The distances at which the ray enters and exits the scene bounding box.
  float tmin = 0, tmax = 0;

  // Use kd tree to find intersection.
  Sphere const *sphere = nullptr;
  if (use_tree) {
    sphere = kdTree.traverse(ray, point, distance_near, distance_far, intersect, tmin, tmax);
  }
  else {
    sphere = brute_force_traverse(ray, point, distance_near, distance_far, intersect, tmin, tmax);
  }

  // If the ray intersected with a sphere, figure out what kind of shading it should have.
  if (sphere) {
    // Check what kind of intersection occured:
    // The intersection occured completely outside the bounding box
    if (distance_far < tmin) {
      return RGB_Dark_Gray;
    }
      // The ray intersects with the sphere before it intersects with the scene bounding box
    else if (chop_spheres && distance_near < tmin) {
      return RGB_Red;
    }
      // The intersection happened within the scene bounding box
    else {
      float norm[3];
      subtractVec(point, sphere->center, norm, 3);
      normalizeVec(norm, 3);
      // Compute color
      float resulting_color[3];
      hadamardVec(sphere->color_reflectivity, light_color, resulting_color, 3);

      // Compute intensity - lambertian shading.
      float intensity = (1. - base_illumination) * clamp(dotVec(norm, light_direction, 3)) + base_illumination;
      // Return the color
      return RGBApixel(
          static_cast<int>(intensity * resulting_color[0] * 255),
          static_cast<int>(intensity * resulting_color[1] * 255),
          static_cast<int>(intensity * resulting_color[2] * 255)
      );
    }
  }
  else {
    return intersect ? RGB_Dark_Gray : RGB_Black;
  }
}

inline Sphere const *RayTrace::brute_force_traverse(const Ray &ray,
                                                    float *point,
                                                    float &distance_near,
                                                    float &distance_far,
                                                    bool &intersect,
                                                    float &tmin,
                                                    float &tmax) const {
  // Find tmin, tmax
  tmin = 0;
  tmax = 0;
  // Check if the ray intersected with the scene bounding box
  intersect = kdTree.getRayIntersectionParameters(ray, tmin, tmax);
  // If the ray does not intersect with the scene bounding box, return.
  if (!intersect) {
    return nullptr;
  }
  // Exhaustively test for ray-sphere intersection
  Sphere const *min_sphere = nullptr;
  float min = 10000, test_point[3];
  for (auto &sphere : sphere_list) {
    if (sphere.intersect(ray, test_point, distance_near, distance_far, tmin) && distance_near < min) {
      min = distance_near;
      min_sphere = &sphere;
      copyVec(test_point, point, 3);
    }
  }
  distance_near = min;
  // This returns nullptr if the ray did not intersect with any sphere
  return min_sphere;
}

inline void RayTrace::createKDTree() {
  // Clear any old kdTree
  if (!kdTree.isCreated()) {
    kdTree.clear();
    // Create the KD Tree
    kdTree.head = new RayKDTreeNode;
    createKDTree_help(0, sphere_list.size() - 1, 0, kdTree.head, kdTree.bounds);
  }
  else {
    kdTree.empty();
  }
  // Insert spheres into the KD tree.
  // MUST use auto& since we want the pointer to the actual sphere.
  for (auto &sphere : sphere_list) {
    kdTree.insert(&sphere);
  }
}

inline void RayTrace::createKDTree_help(int start, int end, int dim, RayKDTreeNode *node, Bounds &bnds) {
  // If there does not need to be any more splitting, return.
  if (end - start + 1 <= target_leaf_size) {
    return;
  }

  // Decide along which dimension to sort the particles - split widest dimension.
  if (adaptive_split) {
    int split_dim = 0;
    if (bnds.wd(0) < bnds.wd(1)) {
      if (bnds.wd(2) < bnds.wd(1)) {
        split_dim = 1;
      }
      else {
        split_dim = 2;
      }
    }
    else {
      if (bnds.wd(2) < bnds.wd(0)) {
        split_dim = 0;
      }
      else {
        split_dim = 2;
      }
    }
    dim = split_dim;
  }
  // If the box is small enough, return.
  if (bnds.wd(dim) <= min_dimension) {
    return;
  }

  // Sort the particles
  quick_sort(start, end, dim);

  // Pick position between center spheres to be the split
  int id = static_cast<int>(0.5 * (start + end));
  float split_value = 0;
  if (id + 1 < sphere_list.size()) {
    split_value = 0.5 * (sphere_list[id].center[dim] + sphere_list[id + 1].center[dim]);
  }
  else {
    return;
  } // This must be the only particle
  // Set split dim and split value
  node->split_dim = dim;
  node->split_value = split_value;

  // Create bounds objects for sub domains
  Bounds lower_bounds(bnds);
  lower_bounds.max[dim] = split_value;
  Bounds upper_bounds(bnds);
  upper_bounds.min[dim] = split_value;

  // Create a left subtree
  node->lower = new RayKDTreeNode;
  createKDTree_help(start, id - 1, (dim + 1) % 3, node->lower, lower_bounds);
  // Create a right subtree
  node->higher = new RayKDTreeNode;
  createKDTree_help(id, end, (dim + 1) % 3, node->higher, upper_bounds);
}

inline void RayTrace::quick_sort(int start, int end, int dim) {
  if (start < end && dim < 3) {
    int partition = quick_sort_partition(start, end, dim);
    quick_sort(start, partition, dim);
    quick_sort(partition + 1, end, dim);
  }
}

inline int RayTrace::quick_sort_partition(int start, int end, int dim) {
  RealType pivot = sphere_list[(start + end) / 2].center[dim];
  int i = start - 1, j = end + 1;
  while (true) {
    do {
      ++i;
    } while (sphere_list[i].center[dim] < pivot);
    do {
      --j;
    } while (sphere_list[j].center[dim] > pivot);
    // Termination condition
    if (j <= i) {
      return j;
    } // Return the id of the new dividing point.
    // If not, swap particles i and j
    swap(sphere_list[i], sphere_list[j]);
  }
}
