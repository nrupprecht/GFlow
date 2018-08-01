#include "clustering.hpp"
// Other files
#include "simdata.hpp"
#include "domainbase.hpp"
#include "printingutility.hpp"

namespace GFlowSimulation {

  Clustering::Clustering(GFlow *gflow) : Base(gflow), same_type_clusters(true), skin(0.) {};

  void Clustering::findClusters() {
    // Look for clusters of particles
    int number = Base::simData->number;
    // If there is nothing to record
    if (number<=0) return;
    // Array used to see which particles are in which cluster - reset it
    // -1 means not on the stack
    clusters = vector<int>(number, -1);
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
      // Go through stack of neighbors, neighbors of neigbors, etc. until none are left
      while (!check_stack.empty()) {
        // ID of the particle to check the type of
        int id = check_stack.top();
        check_stack.pop(); // Pop off the stack
        // Particle with id [id] belongs in the same cluster as particle with id [head]
        clusters[id] = n_clusters;
        // Fill the neighbors array with 
        Base::domain->getAllWithin(id, 2*Base::simData->sg[id] + skin, neighbors);
        for (auto n : neighbors) {
          // Only put on the stack if it is not already on the stack, and has not had its
          if (clusters[n]==-1 && (!same_type_clusters || type[n]==head_type) && type[n]!=-1) {
            check_stack.push(n); // We will look for n's neighbors
            clusters[n] = n_clusters;
          }
        }
      }
      // Increment [head_count]
      ++n_clusters;
      // Get the next head
      get_next(head, clusters);
    }
  }

  vector<int> Clustering::getClusterSizes() const {
    // Do staticstics on cluster size
    vector<int> cluster_size(n_clusters, 0);
    for (int n=0; n<clusters.size(); ++n) ++cluster_size[ clusters[n] ];
    // Return the vector
    return cluster_size;
  }

  const vector<int>& Clustering::getClusters() const {
    return clusters;
  }

  int Clustering::getNumClusters() const {
    return n_clusters;
  }

  void Clustering::get_next(int &head, const vector<int>& cluster) {
    while (head<cluster.size()) {
      if (cluster[head]==-1 && Base::simData->type[head]!=-1) return;
      ++head;
    }
  }

}