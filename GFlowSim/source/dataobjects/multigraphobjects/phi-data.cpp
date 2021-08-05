#include <dataobjects/multigraphobjects/phi-data.hpp>

using namespace GFlowSimulation;

PhiData::PhiData(GFlow *gflow)
    : MultiGraphObject(gflow, "PhiData", "time", "phi", 1) {};

void PhiData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  Vec count_volume(1);
  int size = simData->size_owned();
  for (int n = 0; n < size; ++n) {
    if (-1 < simData->Type(n)) {
      count_volume[0] += sphere_volume(simData->Sg(n), sim_dimensions);
    }
  }

  // Create an entry with the average data. This function handles multiprocessor runs correctly.
  gatherData(gflow->getElapsedTime(), count_volume);

  // If rank 0, correct data so it is a phi, not a total volume. The function getY gets the most recently added datapoint.
  if (topology->getRank() == 0) {
    getY(0) /= (gflow->getBounds().vol() - total_excluded_volume);
  }
}
