#include "clustering.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../base/domainbase.hpp"
#include "../utility/printingutility.hpp"

namespace GFlowSimulation {

  Clustering::Clustering(GFlow *gflow) : Base(gflow), same_type_clusters(true), skin(0.), n_clusters(0), max_cluster_size(0) {};

  void Clustering::findClusters() {
    // Look for clusters of particles
    int number = Base::simData->number;
    // If there is nothing to record
    if (number<=0) return;
    // Array used to see which particles are in which cluster - reset it
    // -1 means not on the stack
    clusters = vector<int>(number, -1);
    cluster_sizes.clear();
    // Get the types
    const int *type = Base::simData->type;
    // Stack for clustering
    std::stack<int> check_stack;
    // Pass into the domain [getAllWithin] function
    vector<int> neighbors;

    // --- Cluster all particles
    int head = 0, head_type(-1);
    // Reset number of clusters
    n_clusters = 0;
    while(head < number) {
      head_type = type[head];
      check_stack.push(head);
      cluster_sizes.push_back(0);
      // Go through stack of neighbors, neighbors of neigbors, etc. until none are left
      while (!check_stack.empty()) {
        // ID of the particle to check the type of
        int id = check_stack.top();
        check_stack.pop(); // Pop off the stack
        ++cluster_sizes[n_clusters]; // Another particle in in the cluster
        // Particle with id [id] belongs in the same cluster as particle with id [head]
        clusters[id] = n_clusters;
        // Fill the neighbors array with ids. The [getAllWithin] function clears [neighbors] before filling it.
        Base::domain->getAllWithin(id, 2*Base::simData->Sg()[id] + skin, neighbors);
        for (auto n : neighbors) {
          // Only put on the stack if it is not already on the stack, and has not had its
          if (clusters[n]==-1 && (!same_type_clusters || type[n]==head_type) && type[n]!=-1) {
            check_stack.push(n); // We will look for n's neighbors
            clusters[n] = n_clusters;
          }
        }
      }
      // Is this cluster the biggest so far?
      if (cluster_sizes.at(n_clusters) > max_cluster_size) 
        max_cluster_size = cluster_sizes.at(n_clusters);
      // Increment [head_count]
      ++n_clusters;
      // Get the next head
      get_next(head, clusters);
    }
  }

  const vector<int>& Clustering::getClusterSizes() const {
    return cluster_sizes;
  }

  const vector<int>& Clustering::getClusters() const {
    return clusters;
  }

  int Clustering::getNumClusters() const {
    return n_clusters;
  }

  int Clustering::getMaxClusterSize() const {
    return max_cluster_size;
  }

  void Clustering::setSkin(RealType s) {
    skin = s;
  }

  void Clustering::setSameTypeClusters(bool s) {
    same_type_clusters = s;
  }

  void Clustering::get_next(int &head, const vector<int>& cluster) {
    while (head<cluster.size()) {
      if (cluster[head]==-1 && Base::simData->type[head]!=-1) return;
      ++head;
    }
  }

}