#include "raytrace.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "palette.hpp" // For color definitions

namespace GFlowSimulation {

  RayTrace::RayTrace() : pix_x(768), pix_y(768), target_leaf_size(50), min_dimension(0.15), use_tree(true), adaptive_split(true) {
    // Camera center
    camera_center[0] = camera_center[1] = camera_center[2] = 0;
    // Camera orientation
    camera_orientation[0] = 1;
    camera_orientation[1] = camera_orientation[2] = 0;
    // Camera right
    camera_up[0] = camera_up[2] = 0; camera_up[1] = 1;
    // Set bitmap size
    image.SetSize(pix_x, pix_y);
    // Create light direction
    light_direction[0] = 1./sqrt(3);
    light_direction[1] = 1./sqrt(3);
    light_direction[2] = 1./sqrt(3);
    normalizeVec(light_direction, 3);
  }

  void RayTrace::initialize() {
    // Make sure camera vectors are good.
    normalizeVec(camera_orientation, 3);
    normalizeVec(camera_up, 3);
    // Correct camera right vector
    float v = dotVec(camera_up, camera_orientation, 3);
    if (v!=0) {
      if (v==1) { // camera_up == camera_orientation
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

    // Timing
    auto start_time = current_time();
    // Pixel for storing color
    RGBApixel color; 
    // Ray orientation, ray direction vectors, accumulator
    float orientation[3], vx[3], vy[3], acc[3]; 
    // Ray trace each pixel
    for (int y=0; y<pix_y; ++y)
      for (int x=0; x<pix_x; ++x) {

        float dx = (float(x) - pix_x/2.f)/float(pix_x);
        float dy = (float(y) - pix_y/2.f)/float(pix_y);

        scalarMultVec(dx, camera_right, vx, 3);
        scalarMultVec(dy, camera_up, vy, 3);

        addVec(camera_orientation, vx, acc, 3);
        addVec(acc, vy, orientation, 3);
        // Normalize the camera vector
        normalizeVec(orientation, 3);
        // Create the ray
        Ray ray(camera_center, orientation);
        // Trace the ray to get the color
        color = trace(ray);
        // Set the color
        image.SetPixel(x, y, color);
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

  void RayTrace::addSphere(const float* center, const float sigma) {
    sphere_list.push_back(Sphere(center, sigma));
  }

  void RayTrace::setCameraPlacement(const float* c) {
    copyVec(c, camera_center, 3);
  }

  void RayTrace::setCameraOrientation(const float* o) {
    copyVec(o, camera_orientation, 3);
  }

  void RayTrace::setBounds(const Bounds& b) {
    kdTree.bounds = b;
  }

  void RayTrace::empty() {
    kdTree.empty();
    sphere_list.clear();
  }

  double RayTrace::getRenderTime() {
    return render_time;
  }

  RGBApixel RayTrace::trace(const Ray& ray) const {
    // Distance to intersection point
    float distance = 0;
    // Intersection point
    float point[3];
    // Did the ray intersect with the scene bounding box?
    bool intersect = false;

    // Use kd tree to find intersection.
    Sphere const* sphere = nullptr;
    if (use_tree) sphere = kdTree.traverse(ray, point, distance, intersect);
    else sphere = brute_force_traverse(ray, point, distance, intersect);

    // If the ray intersected with a sphere, figure out what kind of shading it should have.
    if (sphere) {
      float norm[3];
      subtractVec(point, sphere->center, norm, 3);
      normalizeVec(norm, 3);
      // Compute intensity - lambertian shading.
      float intensity = (1.-base_illumination)*clamp(dotVec(norm, light_direction, 3)) + base_illumination;
      // Return the color
      return RGBApixel(0, intensity*255, 0);
    }
    else return intersect ? RGB_Dark_Gray : RGB_Black;
  }

  inline Sphere const * RayTrace::brute_force_traverse(const Ray& ray, float *point, float& distance, bool& intersect) const {
    Sphere const * min_sphere = nullptr;
    float min = 10000, test_point[3], sphere_center[3];
    bool did_intersect = false;
    for (auto &sphere : sphere_list) {
      if (sphere.intersect(ray, test_point, distance) && distance<min) {
        min = distance;
        min_sphere = &sphere;
        copyVec(test_point, point, 3);
      }
    }
    distance = min;
    // Check if the ray intersected with the scene bounding box
    if (min_sphere) intersect = true;
    else {
      float tmin = 0, tmax = 0;
      intersect = kdTree.getRayIntersectionParameters(ray, tmin, tmax);
    }
    // This returns nullptr if the ray did not intersect with any sphere
    return min_sphere;
  }

  inline void RayTrace::createKDTree() {
    // Clear any old kdTree
    if (!kdTree.isCreated()) {
      kdTree.clear();
      // Create the KD Tree
      kdTree.head = new RayKDTreeNode;
      createKDTree_help(0, sphere_list.size()-1, 0, kdTree.head, kdTree.bounds);
    }
    else kdTree.empty();
    // Insert spheres into the KD tree. 
    // MUST use auto& since we want the pointer to the actual sphere.
    for (auto &sphere : sphere_list) kdTree.insert(&sphere);
  }

  inline void RayTrace::createKDTree_help(int start, int end, int dim, RayKDTreeNode *node, Bounds& bnds) {
    // If there does not need to be any more splitting, return.
    if (end-start+1 <= target_leaf_size) return;
    
    // Decide along which dimension to sort the particles - split widest dimension.
    if (adaptive_split) {
      int split_dim = 0;
      if (bnds.wd(0) < bnds.wd(1)) {
        if (bnds.wd(2)<bnds.wd(1)) split_dim = 1;
        else split_dim = 2;
      }
      else {
        if (bnds.wd(2)<bnds.wd(0)) split_dim = 0;
        else split_dim = 2;
      }
      dim = split_dim;
    }
    // If the box is small enough, return.
    if (bnds.wd(dim) <= min_dimension) return;

    // Sort the particles 
    quick_sort(start, end, dim);
    // Pick position between center spheres to be the split
    int id = static_cast<int>(0.5*(start + end));
    float split_value = sphere_list[id].center[dim];
    node->split_dim = dim;
    node->split_value = split_value;
    
    // Create bounds objects for sub domains
    Bounds lower_bounds(bnds);
    lower_bounds.max[dim] = split_value;
    Bounds upper_bounds(bnds);
    upper_bounds.min[dim] = split_value;

    // Create a left subtree
    node->lower = new RayKDTreeNode;
    createKDTree_help(start, id-1, (dim+1)%3, node->lower, lower_bounds);
    // Create a right subtree
    node->higher = new RayKDTreeNode;
    createKDTree_help(id, end, (dim+1)%3, node->higher, upper_bounds);
  }

  inline void RayTrace::quick_sort(int start, int end, int dim) {
    quick_sort_help(start, end, dim);
  }

  inline void RayTrace::quick_sort_help(int start, int end, int dim) {
    if (start<end && dim<3) {
      int partition = quick_sort_partition(start, end, dim);
      quick_sort_help(start, partition, dim);
      quick_sort_help(partition+1, end, dim);
    }
  }

  inline int RayTrace::quick_sort_partition(int start, int end, int dim) {
    RealType pivot = sphere_list[(start + end)/2].center[dim];
    int i = start-1, j = end+1;
    while (true) {
      do {
        ++i;
      } while (sphere_list[i].center[dim]<pivot);
      do {
        --j;
      } while (sphere_list[j].center[dim]>pivot);      
      // Termination condition
      if (j<=i) return j;
      // If not, swap particles i and j
      swap(sphere_list[i], sphere_list[j]);
    }
  }

  inline void RayTrace::recursion_help(int start, int end, int dim) {
    if (dim>=3) return;
    // Do next level of quicksorts, for the next dimension
    const int sort_bins = 5, min_particles = 10;
    const int ds = (end-start)/sort_bins;

    if (ds>min_particles) {
      for (int i=0; i<sort_bins; ++i) {
        quick_sort_help(i*ds, (i+1)*ds, dim);
        recursion_help (i*ds, (i+1)*ds, dim+1);
      }
      // Potentially, there is some left over
      quick_sort_help(sort_bins*ds, sphere_list.size()-1, dim);
      recursion_help (sort_bins*ds, sphere_list.size()-1, dim+1);
    }
  }

}