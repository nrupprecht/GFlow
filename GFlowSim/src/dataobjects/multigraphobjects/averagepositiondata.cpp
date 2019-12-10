#include "averagepositiondata.hpp"
// Other files
#include "../../base/simdata.hpp"

namespace GFlowSimulation {

  AveragePositionData::AveragePositionData(GFlow *gflow) : MultiGraphObject(gflow, "AvePos", "time", "position", gflow->getSimDimensions()) {
    for (int i=0; i<gflow->getSimDimensions(); ++i)
      axis_y[i] = "Ave pos - X[" + toStr(i) + "]";
  };

  void AveragePositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    Vec ave(sim_dimensions);

    auto x  = simData->X();
    auto im = simData->Im();
    auto type = Base::simData->Type();
    int size = simData->size_owned();
    // Compute totals
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && !isnan(x(n, 0)) && type[n]>-1) { // Presumably, if one component is nan, all are.
        for (int d=0; d<sim_dimensions; ++d)
          ave[d] += x(n, d);
        ++count;
      }

    // Create an entry with the average data. This function handles multiprocessor runs correctly.
    if (useAve) gatherAverageData(Base::gflow->getElapsedTime(), ave, count);
    else gatherData(Base::gflow->getElapsedTime(), ave);
  }

}