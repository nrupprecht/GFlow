#include <dataobjects/positiondata.hpp>
// Other files
#include <utility/printingutility.hpp>
#include <base/simdata.hpp>

using namespace GFlowSimulation;

PositionData::PositionData(GFlow *gflow)
    : DataObject(gflow, "Pos") {
  // The data to gather.
  add_vector_data_entry("X");
  add_magnitude_data_entry("V");
  add_scalar_data_entry("Sg");
  if (1 < gflow->getNTypes()) {
    add_integer_data_entry("Type");
  }
};

void PositionData::pre_integrate() {
  storeData.set_vector_data(vector_data_entries);
  storeData.set_magnitude_data(magnitude_data_entries);
  storeData.set_scalar_data(scalar_data_entries);
  storeData.set_integer_data(integer_data_entries);
  storeData.initialize(simData);
  // Do this after initialize, so bounds are not overwritten.
  storeData.set_data_boundary(gather_bounds);
  // Store initial data (t=0 data).
  storeData.store(initial_data);
}

void PositionData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  // Record and store all the data, optionally using a selection function.
  vector<float> data;
  timeStamps.push_back(gflow->getElapsedTime());
  if (static_cast<bool>(select_function)) {
    storeData.store(data, select_function);
  }
  else {
    storeData.store(data);
  }

  // Collect addition data from other processors.
  MPIObject::mpi_reduce0_position_data(data);
  // If we are the root processor, store the data.
  if (topology->getRank() == 0) {
    positions.push_back(data);
  }
}

bool PositionData::writeToFile(string fileName, bool useName) {
  // The name of the directory for this data
  string dirName = useName ? _correctDirName(fileName) : fileName;
  // Create a directory for all the data
  mkdir(dirName.c_str(), 0777);
  return storeData.write(dirName + dataName + "-" + toStr(object_counter) + ".csv", positions);
}
