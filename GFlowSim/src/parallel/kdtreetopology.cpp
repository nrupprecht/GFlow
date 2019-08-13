#include "kdtreetopology.hpp"

namespace GFlowSimulation {

  KDTreeTopology::KDTreeTopology(GFlow *gflow) : Topology(gflow) {}

  KDTreeTopology::~KDTreeTopology() {
    if (root) delete root;
  }

  //! @brief Compute how the simulation space should be divided up.
  void KDTreeTopology::computeTopology() {
    // Check for valid bounds.
    if (simulation_bounds.vol()<=0) return;

    // Initialize.
    root = new KDTreeTopNode(simulation_bounds, 0);
    compute_decomp(0, numProc, root, 0);

    // Compute neighbors for this processor. 
    find_neighbors();
  }

  //! @brief Given a position and cutoff value, this function returns the
  //! ids of the processors which this particle overlaps.
  void KDTreeTopology::domain_overlaps(const RealType *x, const RealType cutoff, vector<int>& container) {
    // Make sure the container is empty.
    container.clear();

    for (int i=0; i<neighbor_ranks.size(); ++i) {
      if (neighbor_nodes[i]->bounds.distance(x) < cutoff)
        container.push_back(i);
    }
  }

  vector<int> KDTreeTopology::get_neighbor_ranks() const {
    return neighbor_ranks;
  }

  //! @brief Determines which processor a position falls into.
  int KDTreeTopology::domain_ownership(const RealType *x) {
    // First, check if a particle belongs to this domain. As this function is often used to check if particles have left a
    // particular domain, and most will not have, this saves time.
    if (process_bounds.contains(x)) return rank;

    // Otherwise, step through the tree.
    KDTreeTopNode *node = root;
    while (true) {
      // If this is a leaf node.
      if (node->rank!=-1) return node->rank;

      // Otherwise, descend tree.
      if (x[node->split_dim] < node->split_val) node = node->left;
      else node = node->right;
    }
  }

  bool KDTreeTopology::owned_particle(const RealType *x) {
    return process_bounds.contains(x);
  }

  void KDTreeTopology::compute_decomp(int startP, int endP, KDTreeTopNode *node, int dim) {
    // Make sure node is good.
    if (node==nullptr) return;
    // Make sure node has no children.
    if (node->left) {
      delete node->left;
      node->left = nullptr;
    }
    if (node->right) {
      delete node->right;
      node->right = nullptr;
    }

    // Correct dimension
    dim %= sim_dimensions;
    // Assign splitting dimension.
    node->split_dim = dim;

    // Set up
    int nl = static_cast<int>((endP-startP)/2);
    int nr = (endP - startP) - nl;
    Bounds &top_bounds = node->bounds;

    // Divide bounds
    RealType fraction = nl/static_cast<RealType>(endP-startP);
    RealType width = fraction*top_bounds.wd(dim);
    // Set node's splitting value.
    node->split_val = node->bounds.min[dim] + width;

    // Handle left node.
    node->left = new KDTreeTopNode(top_bounds, dim+1 % sim_dimensions);
    node->left->bounds.max[dim] = top_bounds.min[dim] + width; // Adjust bounds
    if (1<nl) compute_decomp(startP, endP-nr, node->left, dim+1);
    else {
      node->left->rank = startP;
      if (startP!=rank) all_processor_nodes.push_back(node->left);
    }
    // Handle right node.
    node->right = new KDTreeTopNode(top_bounds, dim+1 % sim_dimensions);
    node->right->bounds.min[dim] += width;
    if (1<nr) compute_decomp(endP-nr, endP, node->right, dim+1);
    else {
      node->right->rank = endP-1;
      if (endP-1!=rank) all_processor_nodes.push_back(node->right);
    }

    // Check if this process should handle one of these leaves.
    if (node->right && node->right->rank == rank) {
      leaf = node->right;
      process_bounds = node->right->bounds;
    }
    else if (node->left && node->left->rank == rank) {
      leaf = node->left;
      process_bounds = node->left->bounds;
    }

  }

  void KDTreeTopology::find_neighbors() {
    for (auto node : all_processor_nodes) {
      Bounds& bounds = node->bounds;

      // \todo Processors may be neighbors if they are "close enough" to each other, even if they don't directly touch.
      // \todo And what about large particles?

      // Check if the processor bounds are adjacent.
      bool adjacent = true;
      for (int d=0; d<sim_dimensions; ++d) {
        if (bounds.min[d]<=process_bounds.max[d] && process_bounds.min[d]<=bounds.max[d]);
        else {
          adjacent = false;
          break;
        }
      }

      if (adjacent) {
        neighbor_nodes.push_back(node);
        neighbor_ranks.push_back(node->rank);
      }
    }
  }

}