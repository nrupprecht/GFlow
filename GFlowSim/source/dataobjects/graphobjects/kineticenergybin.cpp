#include <dataobjects/graphobjects/kineticenergybin.hpp>

using namespace GFlowSimulation;

KineticEnergyBin::KineticEnergyBin(GFlow *gflow)
    : GraphObject(gflow, "KEBin", "r", "kinetic energy"), center(sim_dimensions) {};

KineticEnergyBin::KineticEnergyBin(GFlow *gflow, RealType r)
    : GraphObject(gflow, "KEBin", "r", "kinetic energy"), radius(r), center(sim_dimensions) {};

void KineticEnergyBin::post_integrate() {
  // Bin width
  RealType dr = radius / bins;
  // Calculate bin distance
  for (int i = 0; i < bins; ++i) {
    data.emplace_back((i + 0.5) * dr, 0.);
  }

  // Bin all kinetic energies.
  auto x = simData->X();
  auto v = simData->V();
  auto im = simData->Im();
  int size = simData->size_owned();
  int sim_dimensions = simData->getSimDimensions();
  for (int n = 0; n < size; ++n) {
    if (im[n] > 0) {
      // Compute distance.
      RealType r = gflow->getDistance(center.data, x(n));
      // Bin data
      if (r <= radius) {
        // Calculate KE
        RealType m = 1. / im[n];
        RealType ke = 0.5 * m * sqr(v(n), sim_dimensions);
        // Find bin
        int b = static_cast<int>(r / dr);
        // Update kinetic energy
        data.at(b).second += ke;
      }
    }
  }
}
