#include "positiondata.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../base/simdata.hpp"
//#include "../visualization/visualization.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow) : DataObject(gflow, "Pos") {
    // The data to gather
    add_vector_data_entry("X");
    add_vector_data_entry("V");
    add_scalar_data_entry("Sg");
    add_integer_data_entry("Type");
    add_scalar_data_entry("StripeX");
    add_scalar_data_entry("Om");
  };

  void PositionData::pre_integrate() {
    storeData.set_vector_data(vector_data_entries);
    storeData.set_scalar_data(scalar_data_entries);
    storeData.set_integer_data(integer_data_entries);
    storeData.initialize(simData);
    // Store initial data
    storeData.store(initial_data);
  }

  void PositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Record what time it was
    float time = Base::gflow->getElapsedTime();
    timeStamps.push_back(time);

    // Record all the data
    vector<float> data;

    // Store the data for this processor.
    storeData.store(data);

    // Collect addition data from other processors.
    MPIObject::mpi_reduce0_position_data(data);
    // If we are the root processor, store the data.
    if (topology->getRank()==0) positions.push_back(data);
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

    return storeData.write(dirName+"data.csv", positions);
  }

  void PositionData::add_vector_data_entry(string entry) {
    vector_data_entries.push_back(entry);
  }

  void PositionData::add_scalar_data_entry(string entry) {
    scalar_data_entries.push_back(entry);
  }

  void PositionData::add_integer_data_entry(string entry) {
    integer_data_entries.push_back(entry);
  }

  void PositionData::clear_all_data_entries() {
    vector_data_entries.clear();
    scalar_data_entries.clear();
    integer_data_entries.clear();
  }

}