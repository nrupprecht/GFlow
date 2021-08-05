#include "dataobjects/multigraphobjects/cavity-statistics.hpp"

using namespace GFlowSimulation;

CavityStatistics::CavityStatistics(GFlow *gflow, real v)
    : MultiGraphObject(gflow, "CavStat", "time", "Number", 3), limit_velocity(v) {
  select_function = [=](const auto& simdata, int n) -> bool {
    return simdata->V(n, 0) < 0.5f * limit_velocity && 0 < simdata->X(n, 0) && 0 < simdata->Im(n);
  };
  storeData.set_vector_data(vector<string>{"X"});
  axis_y[1] = "Average X";
  axis_y[2] = "Average Y";
};

void CavityStatistics::pre_integrate() {
  storeData.initialize(simData);
  // Do this after initialize, so bounds are not overwritten.
  storeData.set_data_boundary(gather_bounds);
}

void CavityStatistics::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  // Record and store all the data, optionally using a selection function.
  vector<float> particle_data;
  storeData.store(particle_data, select_function);
  MPIObject::mpi_reduce0_position_data(particle_data);
  // If we are the root processor, store the data.
  if (topology->getRank() == 0) {
    addEntry(gflow->getElapsedTime());
    int data_width = storeData.getDataWidth(), num_particles = particle_data.size() / data_width;
    getY(0) = num_particles;
    for (int i = 0; i < num_particles; ++i) {
      getY(1) += particle_data[i * data_width];
      getY(2) += particle_data[i * data_width + 1];
    }
    if (0 < num_particles) {
      getY(1) /= num_particles;
      getY(2) /= num_particles;
    }
  }
}
