#include "positiondata.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../base/simdata.hpp"

#include "../visualization/visualization.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow) : DataObject(gflow, "Pos"), dataWidth(0) {
    // The data to gather
    data_types.push_back(DataType::POSITION);
    data_types.push_back(DataType::VELOCITY);
    data_types.push_back(DataType::SIGMA);
    data_types.push_back(DataType::TYPE);
    data_types.push_back(DataType::DISTANCE);
    // Allocate
    vdata = new RealType[sim_dimensions];
  };

  PositionData::~PositionData() {
    if (vdata) delete [] vdata;
  }

  void PositionData::pre_integrate() {
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

  void PositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Record what time it was
    RealType time = Base::gflow->getElapsedTime();
    timeStamps.push_back(time);

    // Record all the data
    vector<RealType> data;
    data.reserve(data_types.size()*simData->Number());
    store_data(data, data_types);
    // Store this timestep's data
    positions.push_back(data);
  }

  bool PositionData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");
    // Make a directory for the data
    mkdir(dirName.c_str(), 0777);
    
    // Print data to csv
    ofstream fout(dirName+"data.csv");
    if (fout.fail()) return false;

    // Print data width, dimensions
    fout << dataWidth << "," << sim_dimensions << "," << positions.size() << "," << Base::simData->ntypes << "\n";

    // Print bounds - mins, then maxes
    Bounds bounds = Base::gflow->getBounds();
    for (int i=0; i<sim_dimensions; ++i) 
      fout << bounds.min[i] << ",";
    for (int i=0; i<sim_dimensions; ++i) {
      fout << bounds.max[i];
      if (i!=sim_dimensions-1) fout << ",";
    }
    fout << "\n";

    // Print out the actual data - first the number of particles, then the particle data
    for (auto &v : positions) fout << v.size() << "," << toCSV(v) << "\n";
    fout.close();

    // Return success
    return true;
  }

  inline void PositionData::store_data(vector<RealType>& data, vector<DataType>& d_types) {
    // Fill the array of data
    for (int i=0; i<Base::simData->Number(); ++i) {
      if (Base::simData->Type(i)!=-1)
        for (const auto type : d_types) get_data(data, type, i);
    }
  }

  inline void PositionData::get_data(vector<RealType>& data, DataType type, int id) {
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
        if (!initial_data.empty()) {
          gflow->getDisplacement(Base::simData->X(id), &initial_data.at(id*sim_dimensions), vdata);        
          data.push_back( magnitudeVec(vdata, sim_dimensions) );
        }
        else data.push_back(-1.);
        break;
      }
    }
  }

}
