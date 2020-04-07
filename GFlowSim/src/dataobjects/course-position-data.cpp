#include "course-position-data.hpp"

namespace GFlowSimulation {

  CoursePositionData::CoursePositionData(GFlow *gflow) : DataObject(gflow, "CPos") {
    // The data to gather.
    add_vector_data_entry("X");
    add_magnitude_data_entry("V");
    add_scalar_data_entry("Sg");
    if (gflow->getNTypes()>1) add_integer_data_entry("Type");
  };

  void CoursePositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Record and store all the data, optionally using a selection function.
    vector<float> data;
    timeStamps.push_back(gflow->getElapsedTime());
    if (static_cast<bool>(select_function)) storeData.store(data, select_function);
    else storeData.store(data);

    // Collect addition data from other processors.
    MPIObject::mpi_reduce0_position_data(data);

    // Discretize data on the root processor.
    if (topology->getRank()==0) {
      
      
    }
  }

}