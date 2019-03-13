#ifndef __RAY_TRACE_HPP__GFLOW__
#define __RAY_TRACE_HPP__GFLOW__

#include "../EasyBMP/EasyBMP.h"
#include "raykdtree.hpp"

namespace GFlowSimulation {

  /*
  struct BSPNode {
    BSPNode() : information(1 << 31), splitCoordinate(0) {};
    
    //unsigned int flagDimAndOffset;
    unsigned int information;
    // If Leaf (flag = 1)
    // bits 0..1      : splitting dimension
    // bits 2..30     : address bits for object list
    // bit  31 (sign) : flag whether node is a leaf
    
    // If Inner (flag = 0)
    // bits 0..1      : empty
    // bits 2..30     : address bits for first son
    // bit  31 (sign) : flat whether node is a leaf
    
    float splitCoordinate;
  };
  inline void setLeafFlag(BSPNode &node, unsigned int flag) {
      node.information = node.information & 0x7FFFFFFF;
      node.information += (flag << 31);
  }
  inline void setAddress(BSPNode &node, unsigned int address) {
      node.information &= 0x80000003;
      node.information += (address << 2);
  }
  inline void setSplitDim(BSPNode &node, unsigned int dim) { node.information += dim; }
  */

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

    //! \brief Place the center of the camera.
    void setCameraPlacement(const float*);

    //! \brief Set the camera's orientation.
    void setCameraOrientation(const float*);

    //! \brief Set the bounds of the KD tree.
    void setBounds(const Bounds&);

    //! \brief Clear the sphere list and KD tree.
    void empty();

  private:
    //! \brief What color a ray should be
    inline RGBApixel trace(const Ray&) const;

    //! \brief Create the KD tree by sorting the spheres.
    inline void createKDTree();

    inline void createKDTree_help(int, int, int, RayKDTreeNode*);

    inline void quick_sort(int, int, int);

    inline void quick_sort_help(int, int, int);

    inline int quick_sort_partition(int, int, int);

    inline void recursion_help(int, int, int);

    //! \brief The target number of spheres to have per leaf node.
    int target_leaf_size;

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
    //! \brief A vector pointint to the "right" relative to the camera.
    float camera_right[3];

    //! \brief The direction of the infinitely distant light source.
    float light_direction[3];

    //! \brief The indirect illumination level. Out of 1.
    float base_illumination = 0.1;

    //! \brief A kd tree for accelerating ray-sphere intersections.
    RayKDTree kdTree;

    //! brief A container for all the spheres in the world.
    vector<Sphere> sphere_list;
  };

}
#endif // __RAY_TRACE_HPP__GFLOW__