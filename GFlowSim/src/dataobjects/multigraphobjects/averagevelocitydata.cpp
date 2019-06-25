#include "averagevelocitydata.hpp"
// Other files
#include "../../base/simdata.hpp"

namespace GFlowSimulation {

  AverageVelocityData::AverageVelocityData(GFlow *gflow) : MultiGraphObject(gflow, "AveVel", "time", "velocity", gflow->getSimDimensions()) {
    for (int i=0; i<gflow->getSimDimensions(); ++i)
      axis_y[i] = "Ave vel - V[" + toStr(i) + "]";
  };

  void AverageVelocityData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    Vec ave(sim_dimensions);

    RealType **v = Base::simData->V();
    RealType *im = Base::simData->Im();
    int size = Base::simData->size(), *type = Base::simData->Type();
    // Compute totals
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && type[n]>-1) {
        for (int d=0; d<sim_dimensions; ++d)
          ave[d] += v[n][d];
        ++count;
      }
    // Check if there was anything to store
    if (count==0) return;

    // Add a new entry to modify
    addEntry(Base::gflow->getElapsedTime());

    // Store data
    RealType c = static_cast<RealType>(count);
    for (int d=0; d<sim_dimensions; ++d) {
      getY(d) = ave[d] / c;
    }
  }

}