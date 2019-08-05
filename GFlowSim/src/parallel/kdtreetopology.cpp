#include "kdtreetopology.hpp"
// Other files
#include "../gflow.hpp"

namespace GFlowSimulation {

  KDTreeTopology::KDTreeTopology(int d) : Topology(d) {}

  KDTreeTopology::~KDTreeTopology() {
    if (root) delete root;
  }

  //! @brief Compute how the simulation space should be divided up.
  void KDTreeTopology::computeTopology() {
    // Check for valid bounds.
    if (simulation_bounds.vol()<=0) return;

    // Initialize.
    int startP = 0, endP = numProc;
    Bounds top_bounds = simulation_bounds;
    int dim = 0;

    root = new KDTreeTopNode(simulation_bounds, 0);
    compute_decomp(startP, endP, root, 0);

    // Compute neighbors for this processor

    // --> FOR NOW: just say everyone is a potential neighbor.
    for (int r=0; r<numProc; ++r)
      if (r!=rank) neighbor_ranks.push_back(r);
  }

  //! @brief Given a position and cutoff value, this function returns the
  //! ids of the processors which this particle overlaps.
  void KDTreeTopology::domain_overlaps(const RealType *x, const RealType sg, vector<int>& container) {
    // Make sure the container is empty().
    container.clear();
    // Go through each dimension.
    for (int d=0; d<sim_dimensions; ++d) {
      if (process_bounds.min[d] + sg < x[d] && process_bounds.max[d] - sg); // The particle in entirely contained within this processor.
      else {
        // \todo Write this function for real. Right now, I'm just adding particle that aren't entirely within their processor to 
        // be ghosts for all other processors.
        for (int i=0; i<neighbor_ranks.size(); ++i) container.push_back(i);
        // Since we already added the particle to be a ghost for all other processors, we can return.
        return;
      }
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

  //! @brief Takes in a processor id and dimension, returns whether there is a domain
  //! "above" it in that dimension.
  bool KDTreeTopology::existsDomainUp(int id, int dim) {
    // STUB
    return false;
  }

  //! @brief Takes in a processor id and dimension, returns whether there is a domain
  //! "below" it in that dimension.
  bool KDTreeTopology::existsDomainDown(int id, int dim) {
    // STUB
    return false;
  }

  //! @brief Get the bounds managed by a processor.
  Bounds KDTreeTopology::getBounds(int rnk) {
    // STUB
    return process_bounds;
  }

  void KDTreeTopology::compute_decomp(int startP, int endP, KDTreeTopNode *node, int dim) {
    // Make sure node is good.
    if (node==nullptr) return;

    // Correct dimension
    dim %= sim_dimensions;

    // Set up
    int nl = static_cast<int>((endP-startP)/2);
    int nr = (endP - startP) - nl;
    Bounds &top_bounds = node->bounds;

    // Divide bounds
    RealType fraction = nl/static_cast<RealType>(endP-startP);
    RealType width = fraction*top_bounds.wd(dim);
    node->split_val = node->bounds.min[dim] + width;

    // Handle left node.
    node->left = new KDTreeTopNode(top_bounds, dim);
    node->left->bounds.max[dim] = top_bounds.min[dim] + width; // Adjust bounds
    if (1<nl) compute_decomp(startP, endP-nr, node->left, dim+1);
    else node->left->rank = startP;
    // Handle right node.
    node->right = new KDTreeTopNode(top_bounds, dim);
    node->right->bounds.min[dim] += width;
    if (1<nr) compute_decomp(endP-nr, endP, node->right, dim+1);
    else node->right->rank = endP-1;

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

}