#include "velocity-volumeplot.hpp"

namespace GFlowSimulation {

  VelocityVolumePlot::VelocityVolumePlot(GFlow *gflow) : VolumePlotObject2D(gflow, "VelocityVolPlot", 2) {
    entry_names[0] = "V[0]";
    entry_names[1] = "V[1]";
  };

  void VelocityVolumePlot::post_step() {
    auto x = simData->X();
    auto v = simData->V();
    Vec V(2);
    for (int i=0; i<simData->size(); ++i) {
      if (focus_bounds.contains(x(i))) {
        // Have the pointer of V point to v.
        V.wrap(v(i), 2);
        // Add to bin by position.
        addToBin(x(i, 0), x(i, 1), V);
        // Stop the pointer of V from pointing to v.
        V.unwrap();
      }
    }
    // Increment
    ++recorded_frames;
  }

}