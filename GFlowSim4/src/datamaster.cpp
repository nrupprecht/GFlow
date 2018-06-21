#include "datamaster.hpp"

namespace GFlowSimulation {

  DataMaster::DataMaster(GFlow *gflow) : Base(gflow) {};

  DataMaster::~DataMaster() {
    for (auto &dob : dataObjects)
      if (dob) delete dob;
  }

  void DataMaster::collect() {
    for (auto dob : dataObjects)
      if (dob) dob->collect();
  }

  void DataMaster::writeToDirectory(string writeDirectory) {
    // --- Do file related things

    // Remove previously existing files if they exist
    system(("rm -rf "+writeDirectory).c_str());

    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);

    // --- Write a summary
    // *** TODO ***

    // --- Have all the data objects write their data
    for (auto dob : dataObjects)
      if (dob) dob->writeToFile(writeDirectory+"/", true);

  }

}