#include "oscillation-data.hpp"

namespace GFlowSimulation {

  OscillationData::OscillationData(GFlow *gflow) : MultiGraphObject(gflow, "Oscillation", "time", "deviation", gflow->getSimDimensions()) {
    for (int i=0; i<gflow->getSimDimensions(); ++i)
      axis_y[i] = "Deviation - X[" + toStr(i) + "]";

    // Hardcode this for now.
    auto bounds = gflow->getBounds();
    for (int d=0; d<sim_dimensions; ++d) {
      real center = 0.5*(bounds.max[d] + bounds.min[d]);
      gather_bounds.min[d] = center - 0.5*0.75*bounds.wd(d);
      gather_bounds.max[d] = center + 0.5*0.75*bounds.wd(d);
    }
  };

  void OscillationData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    Vec pos(sim_dimensions);
    auto x  = simData->X();
    auto im = simData->Im();
    auto type = simData->Type();
    int size = simData->size_owned();

    // Compute totals
    for (int n=0; n<size; ++n) {
      if (im[n]>0 && !isnan(x(n, 0)) && type[n]>-1 && gather_bounds.contains(x(n))) { // Presumably, if one component is nan, all are.
        RealType mass = 1./im[n];
        for (int d=0; d<sim_dimensions; ++d)
          pos[d] += mass*x(n, d);
      }
    }

    // Create an entry with the average data. This function handles multiprocessor runs correctly.
    gatherData(Base::gflow->getElapsedTime(), pos);
  }

  void OscillationData::post_integrate() {
    if (topology->getRank()==0) {
      // Sum.
      Vec ave(sim_dimensions);
      for (int i=0; i<ndata_points; ++i)
      	for (int j=1; j<ndata+1; ++j)
      	  ave[j-1] += multi_data[j][i];
      // Compute average value
      for (int i=0; i<sim_dimensions; ++i) ave[i] /= multi_data.size();
      // Subtract away average value
      for (auto entry : multi_data) {
      	for (int i=1; i<entry.size(); ++i)
      	  entry[i] -= ave[i-1];
      }
    }
  }

}
