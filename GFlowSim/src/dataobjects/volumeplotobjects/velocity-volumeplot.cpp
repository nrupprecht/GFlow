#include "velocity-volumeplot.hpp"

namespace GFlowSimulation {

  VelocityVolumePlot::VelocityVolumePlot(GFlow *gflow) : VolumePlotObject2D(gflow, "VelocityVolPlot", 2) {
    entry_names[0] = "V[0]";
    entry_names[1] = "V[1]";
  };

  void VelocityVolumePlot::post_step() {
    RealType **x = simData->X();
    RealType **v = simData->V();
    Vec V(2);
    for (int i=0; i<simData->size(); ++i) {
      if (focus_bounds.contains(x[i])) {
        V.wrap(v[i], 2);
        addToBin(x[i][0], x[i][1], V);
        V.unwrap();
      }
    }
    // Increment
    ++recorded_frames;
  }

}