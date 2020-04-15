#include "ending-snapshot.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../visualization/visualization.hpp"

namespace GFlowSimulation {

  EndingSnapshot::EndingSnapshot(GFlow *gflow) : DataObject(gflow, "Snapshot") {
    // The data to gather
    vector_data_entries.push_back("X");
    vector_data_entries.push_back("V");
    scalar_data_entries.push_back("Sg");
    scalar_data_entries.push_back("StripeX");
    integer_data_entries.push_back("Type");
    //integer_data_entries.push_back("ID");
  };

  void EndingSnapshot::pre_integrate() {
    storeData.set_vector_data(vector_data_entries);
    storeData.set_scalar_data(scalar_data_entries);
    storeData.set_integer_data(integer_data_entries);
    storeData.initialize(simData);
  }

  void EndingSnapshot::post_integrate() {
    // Gather data from this processor.
    storeData.store(final_data);
    // Collect addition data from other processors.
    MPIObject::mpi_reduce0_position_data(final_data);
    // If we are not the root, erase the data. It is unnecessary.
    if (topology->getRank()!=0) final_data.clear();
  }

  bool EndingSnapshot::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");
    // Make a directory for the data
    mkdir(dirName.c_str(), 0777);
 
    try {
      // Make some images
      Visualization vis;
      vis.setMetaParameters(storeData); // Let the visualization object have the data positions and other meta data
      int dataWidth = storeData.getDataWidth();
      // Kinetic energy snapshot
      vis.setColorOption(2);
      if (sim_dimensions<3) vis.createImage(dirName+"/kinetic.bmp", final_data);
      else vis.createImage3d(dirName+"/kinetic.bmp", final_data);
      // Direction snapshot
      vis.setColorOption(3);
      if (sim_dimensions<3) vis.createImage(dirName+"/direction.bmp", final_data);
      else vis.createImage3d(dirName+"/direction.bmp", final_data);
      // Type snapshot
      if (gflow->getNTypes()>1) {
      	vis.setColorOption(0);
      	if (sim_dimensions<3) vis.createImage(dirName+"/types.bmp", final_data);
      	else vis.createImage3d(dirName+"/types.bmp", final_data);
      }
      // Possible stripex snapshot.
      if (simData->getScalarData("StripeX")!=-1) {
        // Type snapshot
        vis.setColorOption(5);
        if (sim_dimensions<3) vis.createImage(dirName+"/stripe-x.bmp", final_data);
        else vis.createImage3d(dirName+"/stripe-x.bmp", final_data);
      }
    }
    catch (...) {
      cout << "Error in creating snapshot. Continuing.\n";
    }

    // Write the actual data.
    return storeData.write(dirName+"data.csv", final_data);
  }

}
