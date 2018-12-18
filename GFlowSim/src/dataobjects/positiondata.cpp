#include "positiondata.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../base/simdata.hpp"

#include "../visualization/visualization.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow) : DataObject(gflow, "Pos"), dataWidth(0) {
    // The data to gather
    vector_data_entries.push_back("X");
    vector_data_entries.push_back("V");
    scalar_data_entries.push_back("Sg");
    scalar_data_entries.push_back("StripeX");
    integer_data_entries.push_back("Type");
    integer_data_entries.push_back("ID");
    // Allocate
    vdata = new RealType[sim_dimensions];    
  };

  PositionData::~PositionData() {
    if (vdata) delete [] vdata;
  }

  void PositionData::pre_integrate() {
    // Get data positions, calculate data width.
    dataWidth = 0;
    vector<string> temp;
    for (auto entry : vector_data_entries) {
      int pos = simData->get_vector_data(entry);
      if (-1<pos) {
        temp.push_back(entry);
        vector_data_positions.push_back(pos);
        dataWidth += sim_dimensions;
      }
    }
    vector_data_entries = temp;
    temp.clear();
    for (auto entry : scalar_data_entries) {
      int pos = simData->get_scalar_data(entry);
      if (-1<pos) {
        temp.push_back(entry);
        scalar_data_positions.push_back(pos);
        ++dataWidth;
      }
    }
    scalar_data_entries = temp;
    temp.clear();
    for (auto entry : integer_data_entries) {
      int pos = simData->get_integer_data(entry);
      if (-1<pos) {
        temp.push_back(entry);
        integer_data_positions.push_back(pos);
        ++dataWidth;
      }
    }
    integer_data_entries = temp;

    // Store initial data
    initial_data = vector<RealType>(dataWidth*simData->number(), 0);
    store_data(initial_data);
  }

  void PositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Record what time it was
    RealType time = Base::gflow->getElapsedTime();
    timeStamps.push_back(time);

    // Record all the data
    vector<RealType> data(dataWidth*simData->number(), 0);
    store_data(data);
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
    fout << dataWidth << "," << sim_dimensions << "," << positions.size() << "," << Base::simData->ntypes() << "\n";

    // Print bounds - mins, then maxes
    Bounds bounds = Base::gflow->getBounds();
    for (int i=0; i<sim_dimensions; ++i) 
      fout << bounds.min[i] << ",";
    for (int i=0; i<sim_dimensions; ++i) {
      fout << bounds.max[i];
      if (i!=sim_dimensions-1) fout << ",";
    }
    fout << "\n";

    // Vector data types
    fout << vector_data_entries.size() << ",";
    for (int i=0; i<vector_data_entries.size(); ++i) {
      fout << vector_data_entries[i];
      if (i!=vector_data_entries.size()-1) fout << ",";
    }
    fout << "\n";

    // Scalar data types
    fout << scalar_data_entries.size() << ",";
    for (int i=0; i<scalar_data_entries.size(); ++i) {
      fout << scalar_data_entries[i];
      if (i!=scalar_data_entries.size()-1) fout << ",";
    }
    fout << "\n";
    
    // Integer data types
    fout << integer_data_entries.size() << ",";
    for (int i=0; i<integer_data_entries.size(); ++i) {
      fout << integer_data_entries[i];
      if (i!=integer_data_entries.size()-1) fout << ",";
    }
    fout << "\n";

    // Print out the actual data - first the number of particles, then the particle data
    for (auto &v : positions) fout << v.size() << "," << toCSV(v) << "\n";
    fout.close();

    // Return success
    return true;
  }

  inline void PositionData::store_data(vector<RealType>& data) {
    int data_pointer = 0;
    for (int n=0; n<simData->size(); ++n) {
      // If not a particle, continue
      if (simData->Type(n)<0) continue;
      // Copy data
      for (auto v : vector_data_positions) {
        copyVec(simData->VectorData(v)[n], &data[data_pointer], sim_dimensions);
        data_pointer += sim_dimensions;
      }
      for (auto s : scalar_data_positions) {
        data[data_pointer] = simData->ScalarData(s)[n];
        ++data_pointer;
      }
      for (auto i : integer_data_positions) {
        data[data_pointer] = simData->IntegerData(i)[n];
        ++data_pointer;
      }
    }
  }

}
