#include <dataobjects/volumeplotobjects/radius-volumeplot.hpp>

using namespace GFlowSimulation;

RadiusVolumePlot::RadiusVolumePlot(GFlow *gflow)
    : VolumePlotObject2D(gflow, "RadiusVolPlot", 1) {
  entry_names[0] = "R";
};

void RadiusVolumePlot::post_step() {
  auto x = simData->X();
  auto r = simData->Sg();
  for (int i = 0; i < simData->size_owned(); ++i) {
    if (min_radius < r[i] && r[i] < max_radius && gather_bounds.contains(x(i))) {
      addToBin(x(i, 0), x(i, 1), r[i]);
    }
  }
  // Increment
  ++recorded_frames;
}
