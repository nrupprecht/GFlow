#include "averagepositiondata.hpp"
// Other files
#include "../../base/simdata.hpp"

namespace GFlowSimulation {

  AveragePositionData::AveragePositionData(GFlow *gflow) : MultiGraphObject(gflow, "AvePos", "time", "position", 2) {
    // Set counts to not be written
    write_data[1] = false;
  };

  void AveragePositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    Vec ave(sim_dimensions);

    RealType **x = Base::simData->X();
    RealType *im = Base::simData->Im();
    int size = Base::simData->size(), *type = Base::simData->Type();
    // Compute totals
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && type[n]>-1) {
        for (int d=0; d<sim_dimensions; ++d)
          ave[d] += x[n][d];
        ++count;
      }
    // Check if there was anything to store
    if (count==0) return;

    // Add a new entry to modify
    addEntry();

    // Set the time
    getX() = Base::gflow->getElapsedTime();

    // Store data
    for (int d=0; d<sim_dimensions; ++d) {
      ave[d] /= static_cast<RealType>(count);
      getY(d) = ave[d];
    }
  }

}