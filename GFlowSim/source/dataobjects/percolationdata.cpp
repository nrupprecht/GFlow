#include <dataobjects/percolationdata.hpp>

using namespace GFlowSimulation;

PercolationData::PercolationData(GFlow *gflow)
    : DataObject(gflow, "Percolation"), clustering(Clustering(gflow)),
      same_type_clusters(true), skin(0.) {};

PercolationData::PercolationData(GFlow *gflow, RealType x)
    : DataObject(gflow, "Percolation"), clustering(Clustering(gflow)),
      same_type_clusters(true), skin(0.) {};

void PercolationData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  clustering.setSkin(skin);
  clustering.setSameTypeClusters(same_type_clusters);
  clustering.findClusters();
  // Store data
  vector<int> cluster_size = clustering.getClusterSizes();
  cluster_size_record.push_back(cluster_size);
}

bool PercolationData::writeToFile(string fileName, bool useName) {
  // The name of the directory for this data
  string dirName = fileName;
  if (*fileName.rbegin() == '/') { // Make sure there is a /
    dirName += dataName + "/";
  }
  else {
    dirName += ("/" + dataName + "/");
  }

  // Write the data
  // Create a directory for all the data
  mkdir(dirName.c_str(), 0777);
  ofstream fout(dirName + dataName + ".csv");
  if (fout.fail()) {
    return false;
  }
  // Write data here
  for (const auto &clusters : cluster_size_record) {
    fout << toCSV(clusters) << endl;
  }

  // Close file stream
  fout.close();

  // Return success
  return true;
}
