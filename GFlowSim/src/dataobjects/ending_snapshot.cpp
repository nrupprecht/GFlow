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
  };

  void EndingSnapshot::pre_integrate() {
    // Store initial positions
    vector<DataType> pos_data(1, DataType::POSITION);
    store_data(initial_data, pos_data);

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
    final_data.reserve(data_types.size()*simData->Number());
    store_data(final_data, data_types);
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
    BoundsPack boundspack = bounds.pack_up();

    // Make an image
    Visualization vis;
    vis.setColorOption(2);
    vis.createImage(dirName+"/kinetic.bmp", final_data, dataWidth, boundspack, sim_dimensions);
    vis.setColorOption(3);
    vis.createImage(dirName+"/direction.bmp", final_data, dataWidth, boundspack, sim_dimensions);
    vis.setColorOption(4);
    vis.createImage(dirName+"/displacement.bmp", final_data, dataWidth, boundspack, sim_dimensions);
    
    // Print data to csv
    ofstream fout(dirName+"data.csv");
    if (fout.fail()) return false;

    // Print data width, dimensions
    fout << dataWidth << "," << DIMENSIONS << ",1," << Base::simData->ntypes << "\n";

    // Print bounds - mins, then maxes
    for (int i=0; i<DIMENSIONS; ++i) 
      fout << bounds.min[i] << ",";
    for (int i=0; i<DIMENSIONS; ++i) {
      fout << bounds.max[i];
      if (i!=DIMENSIONS-1) fout << ",";
    }
    fout << "\n";

    // Print out the actual data - first the number of particles, then the particle data
    fout << final_data.size() << "," << toCSV(final_data) << "\n";
    fout.close();

    // Return success
    return true;
  }

  inline void EndingSnapshot::store_data(vector<RealType>& data, vector<DataType>& d_types) {
    // Fill the array of data
    for (int i=0; i<Base::simData->Number(); ++i) {
      if (Base::simData->Type(i)!=-1)
        for (const auto type : d_types) get_data(data, type, i);
    }
  }

  inline void EndingSnapshot::get_data(vector<RealType>& data, DataType type, int id) {
    switch(type) {
      case DataType::POSITION: {
        for (int d=0; d<sim_dimensions; ++d) data.push_back(Base::simData->X(id, d));
        break;
      }
      case DataType::VELOCITY: {
        for (int d=0; d<sim_dimensions; ++d) data.push_back(Base::simData->V(id, d));
        break;
      }
      case DataType::SIGMA: {
        data.push_back(Base::simData->Sg(id));
        break;
      }
      case DataType::TYPE: {
        data.push_back(Base::simData->Type(id));
        break;
      }
      case DataType::DISTANCE: {
        RealType displacement[DIMENSIONS];
        gflow->getDisplacement(Base::simData->X(id), &initial_data.at(id*sim_dimensions), displacement);
        data.push_back( magnitudeVec(displacement) );
        break;
      }
    }
  }

}