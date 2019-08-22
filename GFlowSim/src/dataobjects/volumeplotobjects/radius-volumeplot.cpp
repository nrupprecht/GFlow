#include "radius-volumeplot.hpp"

namespace GFlowSimulation {

  RadiusVolumePlot::RadiusVolumePlot(GFlow *gflow) : VolumePlotObject2D(gflow, "RadiusVolPlot", 1) {
    entry_names[0] = "R";
  };

  void RadiusVolumePlot::post_step() {
    RealType **x = simData->X();
    RealType  *r = simData->Sg();
    for (int i=0; i<simData->size(); ++i) {
      if (min_radius < r[i] && r[i] < max_radius && focus_bounds.contains(x[i])) addToBin(x[i][0], x[i][1], r[i]);
    }
    // Increment
    ++recorded_frames;
  }

}