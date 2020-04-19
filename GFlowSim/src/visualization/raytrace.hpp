#ifndef __RAY_TRACE_HPP__GFLOW__
#define __RAY_TRACE_HPP__GFLOW__

#include "../EasyBMP/EasyBMP.h"
#include "raykdtree.hpp"

namespace GFlowSimulation {

  /* 
  *  \brief A class that creates ray-traced 3d images, for system visualization.
  *
  */
  class RayTrace {
  public:
    //! \brief Default constructor.
    RayTrace();

    //! \brief Create necessary structures to do a render.
    void initialize();

    //! \brief Create a bitmap image from the contained data.
    void render();

    //! \brief Save an image to a location.
    void saveImage(std::string);

    //! \brief Reserve some amount of size for the sphere list.
    void reserve(int);

    //! \brief Add a sphere to the world.
    void addSphere(const float*, const float);

    //! \brief Add a sphere to the world, specifying the color.
    void addSphere(const float*, const float, const RGBApixel);

    //! \brief Place the center of the camera.
    void setCameraPlacement(const float*);

    //! \brief Set the camera's orientation.
    void setCameraOrientation(const float*);

    //! \brief Set the bounds of the KD tree.
    void setBounds(const Bounds&);

    //! \brief Set the bounds of the KD tree, scaled so that more (or less) space is included.
    void setBoundsScaled(const Bounds&, real);

    //! \brief Set the resolution of the images the tracer will create.
    void setResolution(int);

    //! \brief Empty the sphere list and KD tree.
    void empty();

    //! \brief Returns the total amount of time the ray tracer has spent rending images.
    double getRenderTime();

  private:
    //! \brief What color a ray should be
    inline RGBApixel trace(const Ray&) const;

    //! \brief Find the sphere that the ray intersects with by brute force.
    //!
    //! Trys intersecting the ray with every sphere in the environment.
    inline Sphere const* brute_force_traverse(const Ray&, float*, float&, float&, bool&, float&, float&) const;

    //! \brief Create the KD tree by sorting the spheres.
    inline void createKDTree();

    //! \brief A helper function for creating the KD tree.
    inline void createKDTree_help(int, int, int, RayKDTreeNode*, Bounds&);

    //! \brief Sort the spheres by the dim-th coordinate using quicksort.
    inline void quick_sort(int, int, int);

    //! \brief Partition into two sections. Return the id of the new dividing point.
    inline int quick_sort_partition(int, int, int);

    //! \brief The target number of spheres to have per leaf node.
    int target_leaf_size;

    //! \brief a box with dimension smaller than this will not be made.
    float min_dimension;

    //! \brief If true, use the kd-tree structure to accelerate ray-sphere intersection. 
    //! Otherwise, use brute force intersection.
    bool use_tree;

    //! \brief If true, the kd tree's split dimension is determined on the fly. If not, the dimensions are
    //! cycled through.
    bool adaptive_split;

    //! \brief If true, chop spheres at the scene bounding box.
    bool chop_spheres;

    //! \brief The bitmap.
    BMP image;

    //! \brief The width of the image.
    int pix_x;
    //! \brief The height of the image.
    int pix_y;

    //! \brief The origin point for the camera.
    float camera_center[3];
    //! \brief The orientation of the camera.
    float camera_orientation[3];
    //! \brief A vector pointing "up" relative to the camera.
    float camera_up[3];

    //! \brief The direction of the infinitely distant light source.
    float light_direction[3];

    //! \brief The color of the infinitely distant light source.
    float light_color[3];

    //! \brief The indirect illumination level. Out of 1.
    float base_illumination = 0.1;

    //! \brief The total amount of time that the tracer has spent rendering.
    double render_time = 0;

    //! \brief A kd tree for accelerating ray-sphere intersections.
    RayKDTree kdTree;

    //! brief A container for all the spheres in the world.
    vector<Sphere> sphere_list;
  };

}
#endif // __RAY_TRACE_HPP__GFLOW__