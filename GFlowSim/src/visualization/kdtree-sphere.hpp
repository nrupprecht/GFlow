#ifndef __KD_TREE_SPHERE_HPP__
#define __KD_TREE_SPHERE_HPP__

#include "../utility/bounds.hpp"

#include "vec3.hpp"

using namespace GFlowSimulation;

struct KDTreeNode {
  //! @brief Destructor, cleans up the subnodes.
  ~KDTreeNode();

  //! @brief The splitting dimension.
  char dim = 0;
  //! @brief The splitting coordinate.
  double split = 0;

  //! @brief The KDTree node for the low side.
  KDTreeNode *min = nullptr;
  //! @brief The KDTree node for the high side.
  KDTreeNode *max = nullptr;

  //! @brief A list of the id's of the spheres in this node. Only for leaf nodes.
  vector<int> contains;
};


class KDTreeSphere {
public:
  ~KDTreeSphere();

  void construct(const vector<Vec3>&, const vector<double>&);

private:
  // --- Helper functions
  inline void constructByHalving();

  KDTreeNode *root = nullptr;

  //! @brief The type of construction the KDTree uses.
  int constructionType = 0;

  BoundsPack bounds;

};

#endif // __KD_TREE_SPHERE_HPP__