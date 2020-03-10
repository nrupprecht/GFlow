#include "average-velocity-data.hpp"
// Other files
#include "../../base/simdata.hpp"

namespace GFlowSimulation {

  AverageVelocityData::AverageVelocityData(GFlow *gflow) : MultiGraphObject(gflow, "AveVel", "time", "velocity", gflow->getSimDimensions()) {
    for (int i=0; i<gflow->getSimDimensions(); ++i)
      axis_y[i] = "Ave vel - V[" + toStr(i) + "]";
    
    // Hardcode this for now.
    auto bounds = gflow->getBounds();
    for (int d=0; d<sim_dimensions; ++d) {
      real center = 0.5*(bounds.max[d] + bounds.min[d]);
      gather_bounds.min[d] = center - 0.5*0.75*bounds.wd(d);
      gather_bounds.max[d] = center + 0.5*0.75*bounds.wd(d);
    }
  };

  void AverageVelocityData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    
    // Get and store data
    Vec ave(sim_dimensions);
    auto x = simData->X();
    auto v = simData->V();
    auto im = simData->Im();
    auto type = simData->Type();
    // Compute totals
    int count = 0;
    for (int n=0; n<simData->size_owned(); ++n)
      if (im[n]>0 && type[n]>-1) {
        if (!isnan(v(n, 0)) && gather_bounds.contains(x(n))) { // Presumably, if one component is nan, all are.
          for (int d=0; d<sim_dimensions; ++d)
            ave[d] += v(n, d);
          ++count;
        }
      }

    // Create an entry with the average data. This function handles multiprocessor runs correctly.
    gatherAverageData(Base::gflow->getElapsedTime(), ave, count);
  }

}