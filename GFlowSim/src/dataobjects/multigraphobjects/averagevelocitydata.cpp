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

    auto v = simData->V();
    auto im = simData->Im();
    auto type = simData->Type();
    int size = simData->size();
    // Compute totals
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && type[n]>-1) {
        if (!isnan(v(n, 0))) { // Presumably, if one component is nan, all are.
          for (int d=0; d<sim_dimensions; ++d)
            ave[d] += v(n, d);
          ++count;
        }
      }
    // Check if there was anything to store
    if (count==0) return;

    // Create an entry with the average data. This function handles multiprocessor runs correctly.
    gatherAverageData(Base::gflow->getElapsedTime(), ave, count);
  }

}