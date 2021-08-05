#include <dataobjects/graphobjects/average-speed-data.hpp>
// Other files
#include <base/simdata.hpp>
#include <utility/vectormath.hpp>
#include <visualization/palette.hpp>

using namespace GFlowSimulation;

AverageSpeedData::AverageSpeedData(GFlow *gflow)
    : GraphObject(gflow, "AveV", "time", "average velocity") {};

void AverageSpeedData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  // Get and store data
  RealType av = 0;
  auto v = simData->V();
  auto im = simData->Im();
  auto type = simData->Type();
  int count = 0;
  for (int n = 0; n < simData->size_owned(); ++n) {
    if (im[n] > 0 && type[n] > -1) {
      av += magnitudeVec(v(n), sim_dimensions);
      ++count;
    }
  }

  // Store data. These functions work correctly with multiprocessor runs.
  gatherAverageData(gflow->getElapsedTime(), av, count);
}
