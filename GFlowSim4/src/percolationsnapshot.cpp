#include "percolationsnapshot.hpp"

namespace GFlowSimulation {

  PercolationSnapshot::PercolationSnapshot(GFlow *gflow) : DataObject(gflow, "PercolationSnapshot"), 
    clustering(Clustering(gflow)) {};

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
      constituents[clus].push_back(Base::simData->sg[i]);
    }

    // Fill [record] vector
    vector<RealType*> record(n_clusters, nullptr);
    for (int i=0; i<n_clusters; ++i) {
      record[i] = new RealType[ constituents[i].size() ];
      copyVec(&constituents[i][0], record.at(i), constituents[i].size());
    }
    // Fill the [elements] vector
    for (int i=0; i<n_clusters; ++i) 
      elements[i] = constituents[i].size();

  }
  


  bool PercolationSnapshot::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    if (!PrintingUtility::writeVectorToDirectory(record, elements, DIMENSIONS + 1, dirName, dataName)) return false;

    return true;
  }

  void PercolationSnapshot::clearRecord() {
    for (auto &ptr : record) {
      if (ptr) delete [] ptr;
    }
  }

}