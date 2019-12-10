#include "avevelocitydata.hpp"
// Other files
#include "../../base/simdata.hpp"
#include "../../utility/vectormath.hpp"
#include "../../visualization/palette.hpp"

namespace GFlowSimulation {
  // Constructor
  AveVelocityData::AveVelocityData(GFlow *gflow) : GraphObject(gflow, "AveV", "time", "average velocity") {};

  void AveVelocityData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType av = 0;
    auto v = Base::simData->V();
    auto im = Base::simData->Im();
    auto type = Base::simData->Type();
    int size = Base::simData->size();
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && type[n]>-1) {
        av += magnitudeVec(v(n), sim_dimensions);
        ++count;
      }

    // Store data. These functions work correctly with multiprocessor runs.
    gatherAverageData(gflow->getElapsedTime(), av, count);
  }

}