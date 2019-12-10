#include "windtunnel.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  WindTunnel::WindTunnel(GFlow *gflow, RealType v) : Modifier(gflow), velocity(v), halfWidth(4.), acceleration(2.) {
    leftBound  = gflow->getBounds().min[0] + halfWidth;
    rightBound = gflow->getBounds().max[0] - halfWidth;
  }

  void WindTunnel::post_forces() {
    // Only run the wind tunnel during an actual simulation.
    if (gflow->getRunMode()!=RunMode::SIM) return;

    // Currently hardcoded to only be able to handle 4 dimensions or fewer.
    if (sim_dimensions>4) throw BadDimension();
        
    auto x = simData->X(), v = simData->V(), f = simData->F();
    int size = simData->size();
    for (int i=0; i<size; ++i) {
      if (x(i, 0) < leftBound || rightBound < x(i,0)) {
        // Act like an overdamped integrator. Hopefully, this will reduce the occurence of waves getting propagated around the tube.
        zeroVec(v(i), sim_dimensions);
        v(i, 0) = velocity;
        plusEqVecScaled(v(i), f(i), 0.25f*DEFAULT_DAMPING_CONSTANT, sim_dimensions);
        // Then, zero the force.
        zeroVec(f(i), sim_dimensions);
      }
    }
  }

}