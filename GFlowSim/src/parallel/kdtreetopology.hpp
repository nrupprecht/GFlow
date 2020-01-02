#ifndef __KD_TREE_TOPOLOGY_HPP__GFLOW__
#define __KD_TREE_TOPOLOGY_HPP__GFLOW__

#include "../base/topology.hpp"

#include "../base/simdata.hpp"

namespace GFlowSimulation {

  struct KDTreeTopNode {
    //! \brief Default constructor, sets bounds.
    KDTreeTopNode(Bounds &b) : bounds(b), split_val(0) {};

    //! \brief Constructor that sets bounds and splitting dimension.
    KDTreeTopNode(Bounds &b, int sd) : bounds(b), split_dim(sd) {};

    //! \brief Destructor.
    ~KDTreeTopNode() {
      if (left) delete left;
      if (right) delete right;
    }

    void print(std::ostream& out, int indent=0) {
      if (rank!=-1) out << string(indent, ' ') << rank << ": " << bounds << "\n";
      else {
        out << string(indent, ' ') << "Node: " << split_dim << ", " << split_val << " :: " << bounds << endl;
        if (left)   left->print(out, indent+1);
        if (right) right->print(out, indent+1);
      }
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


  /* Check if the domains are either the split along dimensions d, or within cutoff distance.
    Spit.
       B|
    ----|----
        |A

    Within cutoff (as long as the distance from A to B is within the cutoff).

       B|
    ----|
        |-----
        |A
  */

  class KDTreeTopology : public Topology {
  public:
    //! \brief Default constructor, takes the number of dimensions.
    KDTreeTopology(GFlow*);

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

    //! \brief Return the bounds of the i-th neighboring processor.
    virtual const Bounds& get_neighbor_bounds(int) const override;

    // --- Object exchange

    //! \brief Exchange particles that belong to other processors.
    //!
    //! Move particles that belong to other domains to those domains, and delete them from here. Then recieve
    //! particles from other domains that belong to this domain.
    virtual void exchange_particles() override;

    //! \brief Send and receive ghost particles.
    virtual void create_ghost_particles() override;

    //! \brief Send data to other processors to update ghost data.
    virtual void send_ghost_updates() override;

    //! \brief Receive data from other processors to update ghosts on this processor.
    virtual void recv_ghost_updates() override;

  private:

    //! \brief Compute the KD tree decomposition of the simulation bounds.
    void compute_decomp(int, int, KDTreeTopNode*, int);

    //! \brief Find the neighbors of this processor.
    //!
    //! This should be called after compute_decomp.
    void find_neighbors();

    //! \brief Determine which neighbors we might have to send ghosts particles to.
    void determine_sendable_neighbor_ranks();

    // --- Sending helpers

    //! \brief Pack up the particle data for the specified ids and send it to another processor, optionally deleting the original particles
    //! from this processor.
    template<unsigned=0> void send_particle_data(const vector<int>&, int, vector<RealType>&, MPI_Request*, int, bool=false);

    //! \brief Send particles, so that their position relative to the other processor is minimal. Used for sending ghost particles.
    template<unsigned=0> void send_particle_data_relative(const vector<int>&, int, vector<RealType>&, MPI_Request*, int, int);

    //! \brief Recieve particle information, and use it to create new particles. Return the number of particles recieved.
    template<unsigned=0> int recv_new_particle_data(int, vector<RealType>&, int);

    // --- Data

    //! \brief The root of the kd tree.
    KDTreeTopNode *root = nullptr;
    //! \brief The node in the tree corresponding to the bounds for this processor.
    KDTreeTopNode *leaf = nullptr;

    //! \brief A vector of all the (other) processor nodes (leaves of the tree).
    vector<KDTreeTopNode*> all_processor_nodes;

    //! \brief The nodes of the i-th neighbors for this node. Useful for checking what neighbor's bounds are.
    vector<KDTreeTopNode*> neighbor_nodes;
  };

  // Include template functions.
  #include "kdtreetopology.tpp"

}

#endif // __KD_TREE_TOPOLOGY_HPP__GFLOW__