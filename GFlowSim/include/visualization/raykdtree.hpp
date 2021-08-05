#ifndef __RAY_KD_TREE_HPP__GFLOW__
#define __RAY_KD_TREE_HPP__GFLOW__

#include "raytracestructs.hpp"

namespace GFlowSimulation {

  struct RayKDTreeNode {
    //! \brief Clean up this nodes's children.
    ~RayKDTreeNode();

    //! \brief Insert a sphere into the KD tree.
    void insert(Sphere*);

    //! \brief Empty the node and all its children of spheres, but keep the same structure.
    void empty();

    //! \brief Traverse the tree, searching for the first sphere the ray intersects with
    Sphere* traverse(const Ray&, float*, float&, float&, float, float) const;

    // --- Accessor helpers 

    //! \brief Get the maximum depth of the KD tree.
    int getMaxDepth() const;

    //! \brief Get the total number spheres on all the nodes, including duplicates.
    int getTotalSpheres() const;

    // --- Data

    //! \brief Along what dimension is the splitting for this tree's children.
    int split_dim = 0;

    //! \brief The value at which this subvolume splits.
    float split_value = 0;

    //! \brief Whether or not this is a terminal node.
    bool terminal = true;

    //! \brief The children of the KD Tree.
    RayKDTreeNode *lower = nullptr, *higher = nullptr;

    //! \brief A vector of pointers to spheres in this KDTree node.
    vector<Sphere*> sphere_list;
    
  };

  class RayKDTree {
  public:
    //! \brief Clean up the tree.
    ~RayKDTree();

    //! \brief Insert a sphere into the KD tree.
    void insert(Sphere*);

    //! \brief Send a ray through the data structure, return the first sphere hit by the Ray, 
    //! and nullptr if no sphere is encountered.
    //!
    //! The boolean is set to whether the ray intersected with the scene bounding box.
    //!
    //! \param ray The ray that is marched through the tree.
    //! \param point The intersection point of the ray with the first sphere it hits.
    //! \param distance_near The distance at which the ray enters the sphere.
    //! \param distance_far The distance at which the ray exits the sphere.
    //! \param intersect Whether or not the ray intersects with the scene bounding box.
    //! \param tmin The distance at which the ray enters the bounding box (0 if the ray's origin is within the bounding box).
    //! \param tmax The distance at which the ray leaves the bounding box (0 if the ray never intersects with the bounding box).
    Sphere* traverse(const Ray&, float*, float&, float&, bool&, float&, float&) const;

    //! \brief Get the part of the ray that intersects with the scene bounding box.
    //!
    //! \param ray The ray in question.
    //! \param tmin The distance at which the ray enters the bounding box (0 if the ray's origin is within the bounding box).
    //! \param tmax The distance at which the ray leaves the bounding box (0 if the ray never intersects with the bounding box).
    bool getRayIntersectionParameters(const Ray&, float&, float&) const;

    //! \brief Clean up the tree.
    void clear();

    //! \brief Empty the tree of spheres, but keep the same structure.
    void empty();

    //! \brief Returns true if the tree has been created.
    bool isCreated();

    // --- Accessors

    //! \brief Get the maximum depth of the KD tree.
    int getMaxDepth() const;

    //! \brief Get the total number spheres on all the nodes, including duplicates.
    int getTotalSpheres() const;

    // Ray Trace is a friend class.
    friend class RayTrace;

  private:
    //! \brief The head of the KD Tree.
    RayKDTreeNode *head = nullptr;

    //! \brief The bounds of the region enclosed in the KD tree.
    Bounds bounds = Bounds(3);
  };

}
#endif // __RAY_KD_TREE_HPP__GFLOW__