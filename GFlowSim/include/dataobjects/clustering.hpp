#ifndef __CLUSTERING_HPP__GFLOW__
#define __CLUSTERING_HPP__GFLOW__

#include "../gflow.hpp"
// For neighbor traversal
#include <stack>

namespace GFlowSimulation {

  class Clustering : public Base {
  public:
    //! Constructor
    Clustering(GFlow*);

    //! Find clusters of particles
    void findClusters();

    //! Get a vector of all the cluster sizes
    const vector<int>& getClusterSizes() const;

    //! Get the clusters vector
    const vector<int>& getClusters() const;

    //! Get the number of clusters
    int getNumClusters() const;

    int getMaxClusterSize() const;

    // --- Mutators

    void setSkin(RealType);

    void setSameTypeClusters(bool);

  private:
    // --- Helper functions

    //! Get the next uncatagorized particle id
    void get_next(int&, const vector<int>&);

    // --- Data

    //! A mapping from particle id to cluster #
    vector<int> clusters;
    //! The sizes of the n-th cluster
    vector<int> cluster_sizes;
    //! The number of clusters
    int n_clusters, max_cluster_size;

    //! If true, we cluster particles of the same type. Otherwise, cluster
    //! all types of particles together
    bool same_type_clusters;

    //! How close to touching should particles be to be counted as being in the same cluster
    RealType skin;
  };
}
#endif // __CLUSTERING_HPP__GFLOW__