#include "percolationsnapshot.hpp"

namespace GFlowSimulation {

  PercolationSnapshot::PercolationSnapshot(GFlow *gflow) : DataObject(gflow, "PercolationSnapshot"), skin(0), same_type_clusters(true), 
    clustering(Clustering(gflow)) {
      clustering.setSkin(skin);
      clustering.setSameTypeClusters(same_type_clusters);
  };

  PercolationSnapshot::PercolationSnapshot(GFlow *gflow, RealType s) : DataObject(gflow, "PercolationSnapshot"), skin(s), same_type_clusters(true), 
    clustering(Clustering(gflow)) {
      clustering.setSkin(skin);
      clustering.setSameTypeClusters(same_type_clusters);
  };

  PercolationSnapshot::~PercolationSnapshot() {
    clearRecord();
  }

  void PercolationSnapshot::post_integrate() {
    // Take a snapshot of the percolation
    clustering.findClusters();

    // Get the particle clustering data
    const vector<int>& clusters = clustering.getClusters();
    // How many clusters are there
    int n_clusters = clustering.getNumClusters();

    vector<RealType> *constituents = new vector<RealType> [ n_clusters ];
    for (int i=0; i<clusters.size(); ++i) {
      int clus = clusters[i];
      // Push the position
      for (int d=0; d<DIMENSIONS; ++d)
        constituents[clus].push_back(Base::simData->X(i,d));
      // Push the radius
      constituents[clus].push_back(Base::simData->Sg()[i]);
    }

    // Fill [record] and [elements] vectors
    clearRecord();
    elements.clear();
    record = vector<RealType*>(n_clusters, nullptr);
    for (int i=0; i<n_clusters; ++i) {
      record[i] = new RealType[ constituents[i].size() ];
      copyVec(&constituents[i][0], record.at(i), constituents[i].size());
      elements.push_back(constituents[i].size()/(DIMENSIONS+1));
    }

    // Clean up
    delete [] constituents;
  }
  


  bool PercolationSnapshot::writeToFile(string fileName, bool useName) {
    if (!PrintingUtility::writeVectorToDirectory(record, elements, DIMENSIONS + 1, fileName, dataName)) return false;
    return true;
  }

  void PercolationSnapshot::clearRecord() {
    for (auto &ptr : record) {
      if (ptr) delete [] ptr;
    }
  }

}