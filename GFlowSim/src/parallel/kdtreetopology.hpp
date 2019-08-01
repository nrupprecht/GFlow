#ifndef __KD_TREE_TOPOLOGY_HPP__GFLOW__
#define __KD_TREE_TOPOLOGY_HPP__GFLOW__

#include "topology.hpp"

namespace GFlowSimulation {

  struct KDTreeTopNode {
    //! \brief Constructor.
    KDTreeTopNode(Bounds &b, int sd) : bounds(b), split_dim(sd) {};

    //! \brief Destructor.
    ~KDTreeTopNode() {
      if (left) delete left;
      if (right) delete right;
    }

    //! \brief What axis we are splitting along.
    int split_dim = 0;
    //! \brief The split value.
    RealType split_val = 0;

    // \brief The rank of the processor that this node corresponds to, if this is a leaf node. If this is not a leaf node, this is set to -1.
    int rank = -1;
    //! \brief The bounds contained in this node.
    Bounds bounds;
    //! \brief Pointers
    KDTreeTopNode *left = nullptr, *right = nullptr;
  };

  class KDTreeTopology : public Topology {
  public:
    //! \brief Default constructor, takes the number of dimensions.
    KDTreeTopology(int);

    //! \brief Destructor.
    virtual ~KDTreeTopology();

    //! \brief Compute how the simulation space should be divided up.
    virtual void computeTopology() override;

    //! \brief Given a position and cutoff value, this function returns the 
    //! ids of the processors which this particle overlaps.
    virtual void domain_overlaps(const RealType*, const RealType, vector<int>&) override;

    //! \brief Return a vector of the ranks of processors that are potential neighbors for this processor.
    virtual vector<int> get_neighbor_ranks() const override;

    //! \brief Determines which processor a position falls into.
    virtual int domain_ownership(const RealType*) override;

    //! \brief Whether this particle should be owned by this processor.
    virtual bool owned_particle(const RealType*) override;

    //! \brief Takes in a processor id and dimension, returns whether there is a domain
    //! "above" it in that dimension.
    virtual bool existsDomainUp(int, int) override;

    //! \brief Takes in a processor id and dimension, returns whether there is a domain
    //! "below" it in that dimension.
    virtual bool existsDomainDown(int, int) override;

    //! \brief Get the bounds managed by a processor.
    Bounds getBounds(int) override;

  private:

    void compute_decomp(int, int, KDTreeTopNode*, int);

    //! \brief The root of the kd tree.
    KDTreeTopNode *root = nullptr;
    //! \brief The node in the tree corresponding to the bounds for this processor.
    KDTreeTopNode *leaf = nullptr;

    //! \brief The ranks of processors that are potential neighbors of this processor.
    vector<int> neighbor_ranks;
  };

}

#endif // __KD_TREE_TOPOLOGY_HPP__GFLOW__