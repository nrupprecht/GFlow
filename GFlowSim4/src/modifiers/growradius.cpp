#include "growradius.hpp"

namespace GFlowSimulation {

  GrowRadius::GrowRadius(GFlow *gflow) : Modifier(gflow) {
    time0 = gflow->getElapsedTime();
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
    simData->Sg(id) = sigma;
  }

}