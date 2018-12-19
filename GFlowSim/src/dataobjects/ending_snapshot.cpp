#include "ending_snapshot.hpp"
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
    integer_data_entries.push_back("ID");
  };

  void EndingSnapshot::pre_integrate() {
    storeData.set_vector_data(vector_data_entries);
    storeData.set_scalar_data(scalar_data_entries);
    storeData.set_integer_data(integer_data_entries);
    storeData.initialize(simData);
  }

  void EndingSnapshot::post_integrate() {
    storeData.store(final_data);
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

    // Get bounds
    Bounds bounds = Base::gflow->getBounds();

    // Make an image
    Visualization vis;
    vis.setMetaParameters(storeData); // Let the visualization object have the data positions and other meta data
    int dataWidth = storeData.getDataWidth();
    vis.setColorOption(2);
    vis.createImage(dirName+"/kinetic.bmp", final_data);
    vis.setColorOption(3);
    vis.createImage(dirName+"/direction.bmp", final_data);
    
    return storeData.write(dirName+"data.csv", final_data);
  }

}