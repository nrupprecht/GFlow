#include "datamaster.hpp"

namespace GFlowSimulation {

  DataMaster::DataMaster(GFlow *gflow) : Base(gflow) {};

  DataMaster::~DataMaster() {
    for (auto &dob : dataObjects)
      if (dob) delete dob;
  }

  void DataMaster::pre_step() {
    for (auto dob : dataObjects)
      if (dob) dob->pre_step();
  }

  void DataMaster::addDataObject(DataObject *dob) {
    dataObjects.push_back(dob);
  }
  
  void DataMaster::pre_exchange() {
    for (auto dob : dataObjects)
      if (dob) dob->pre_exchange();
  }
  
  void DataMaster::pre_neighbors() {
    for (auto dob : dataObjects)
      if (dob) dob->pre_neighbors();
  }

  void DataMaster::pre_forces() {
    for (auto dob : dataObjects)
      if (dob) dob->pre_forces();
  }

  void DataMaster::post_forces() {
    for (auto dob : dataObjects)
      if (dob) dob->post_forces();
  }

  void DataMaster::post_step() {
    for (auto dob : dataObjects)
      if (dob) dob->post_step();
  }

  bool DataMaster::writeToDirectory(string writeDirectory) {
    // --- Do file related things

    // Remove previously existing files if they exist
    system(("rm -rf "+writeDirectory).c_str());

    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);
    // Create a subdirectory for stat data
    mkdir((writeDirectory+"/StatData").c_str(), 0777);

    // --- Write a summary
    // *** TODO ***

    // --- Have all the data objects write their data
    bool success = true;
    for (auto dob : dataObjects)
      if (dob) success &= dob->writeToFile(writeDirectory+"/StatData/", true);

    // Return true if all writes were successful
    return success;
  }

}