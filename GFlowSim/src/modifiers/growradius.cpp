#include "growradius.hpp"

namespace GFlowSimulation {

  GrowRadius::GrowRadius(GFlow *gflow, int global, RealType r0, RealType rf, RealType time) : Modifier(gflow), global_id(global) {
    if (r0<0) throw BadRadius();
    time0 = gflow->getElapsedTime();
    sigma0 = r0;
    sigmaf = rf;
    rdot = (rf-r0)/time;
  };

  void GrowRadius::post_forces() {
    RealType time = gflow->getElapsedTime();
    // Find the index of the particle
    int id = simData->getLocalID(global_id);
    // Calculate what the radius should be
    RealType sigma = (time - time0)*rdot + sigma0;
    if (sigma>sigmaf) {
      sigma = sigmaf;
      // This modifier is over
      remove = true;
    }
    // Set the radius
    if (id>-1) simData->Sg()[id] = sigma;
    else remove = true;
  }

}