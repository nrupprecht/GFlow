#ifndef __KD_TREE_SPHERE_HPP__
#define __KD_TREE_SPHERE_HPP__

struct KDTreeNode {
  //! @brief The splitting dimension.
  char dim;
  //! @brief The splitting coordinate.
  double split;

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

  KDTreeNode *root = nullptr;

};

#endif // __KD_TREE_SPHERE_HPP__