#include "stripe-x.hpp"
// Other files
#include "../base/topology.hpp"

namespace GFlowSimulation {

  StripeX::StripeX(GFlow *gflow) : Modifier(gflow), window(0.5), lastUpdate(0), updateDelay(0.1) {};

  void StripeX::pre_integrate() {
    // Get the strip x entry.
    entry = simData->requestScalarData("StripeX");
    // Entry must be positive. Otherwise, something is wrong.
    if (entry<0) throw false;
    // Initialize all values
    auto st = simData->ScalarData(entry);
    auto x = simData->X();
    // Set initial heights
    for (int i=0; i<simData->size_owned(); i+=sim_dimensions)
      st(i) = x(i, 1);
  }
  
  void StripeX::post_forces() {

    return;

    // Get the time
    RealType time = gflow->getElapsedTime();
    if (time-lastUpdate<updateDelay || entry<0) return;

    // Set particles in the window. We only need to do something if this processor overlaps with these bounds.
    RealType bound = simData->getBounds().min[0] + window;
    if (bound<topology->getProcessBounds().min[0]) return;

    // Set heights for particles in the window
    auto x = simData->X();
    auto st = simData->ScalarData(entry);
    for (int i=0; i<simData->size_owned(); i+=sim_dimensions)
      if (x(i, 0) < bound) st(i) = x(i, 1);

    // Update time point
    lastUpdate = time;
  }

}