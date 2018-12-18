#include "ending_snapshot.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../visualization/visualization.hpp"

namespace GFlowSimulation {
  // Constructor
  EndingSnapshot::EndingSnapshot(GFlow *gflow) : DataObject(gflow, "Snapshot"), dataWidth(0) {
    // The data to gather
    data_types.push_back(DataType::POSITION);
    data_types.push_back(DataType::VELOCITY);
    data_types.push_back(DataType::SIGMA);
    data_types.push_back(DataType::TYPE);
    data_types.push_back(DataType::DISTANCE);
    // Allocate
    vdata = new RealType[sim_dimensions];
  };

  EndingSnapshot::~EndingSnapshot() {
    if (vdata) delete [] vdata;
  }

  void EndingSnapshot::pre_integrate() {
    // Store initial positions
    vector<DataType> pos_data(1, DataType::POSITION);
    store_data<RealType>(initial_data, pos_data);

    // Calculate data width
    dataWidth = 0;
    for (const auto type : data_types) {
      if (type==DataType::POSITION || type==DataType::VELOCITY)
        dataWidth += sim_dimensions;
      else ++dataWidth;
    }
  }

  void EndingSnapshot::post_integrate() {
    // Record all the data
    final_data.reserve(data_types.size()*simData->number());
    store_data<double>(final_data, data_types);
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
    vis.setColorOption(2);
    vis.createImage(dirName+"/kinetic.bmp", final_data, dataWidth, bounds, sim_dimensions);
    vis.setColorOption(3);
    vis.createImage(dirName+"/direction.bmp", final_data, dataWidth, bounds, sim_dimensions);
    vis.setColorOption(4);
    vis.createImage(dirName+"/displacement.bmp", final_data, dataWidth, bounds, sim_dimensions);
    
    // Print data to csv
    ofstream fout(dirName+"data.csv");
    if (fout.fail()) return false;

    // Print data width, dimensions
    fout << dataWidth << "," << sim_dimensions << ",1," << Base::simData->ntypes() << "\n";

    // Print bounds - mins, then maxes
    for (int i=0; i<sim_dimensions; ++i) 
      fout << bounds.min[i] << ",";
    for (int i=0; i<sim_dimensions; ++i) {
      fout << bounds.max[i];
      if (i!=sim_dimensions-1) fout << ",";
    }
    fout << "\n";

    // Print out the actual data - first the number of particles, then the particle data
    fout << final_data.size() << "," << toCSV(final_data) << "\n";
    fout.close();

    // Return success
    return true;
  }

}