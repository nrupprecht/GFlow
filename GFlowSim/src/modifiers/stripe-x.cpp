#include "stripe-x.hpp"
// Other files
#include "../base/interactionhandler.hpp"
#include "../base/topology.hpp"

namespace GFlowSimulation {

  StripeX::StripeX(GFlow *gflow) : Modifier(gflow), window(0.5), lastUpdate(0), updateDelay(0.1) {};

  void StripeX::pre_integrate() {
    if (sim_dimensions==1) return;
    // Get the strip x entry.
    entry = simData->requestScalarData("StripeX");

    // Entry must be positive. Otherwise, something is wrong.
    if (entry<0) throw false;
    // Initialize all values
    auto st = simData->ScalarData(entry);
    auto x = simData->X();

    // Set initial heights
    for (int i=0; i<simData->size_owned(); ++i)
      st[i] = x[i][1];
  }
  
  void StripeX::post_forces() {
    if (sim_dimensions==1) return;
    // Get the time
    RealType time = Base::gflow->getElapsedTime();
    if (time-lastUpdate<updateDelay || entry<0) return;
    // Update time point
    lastUpdate = time;

    // Set particles in the window
    RealType bound = handler->getSimulationBounds().min[0] + window;

    // Check if this processor can have any particles below the bound.
    if (bound<handler->getProcessBounds().min[0]) return;

    auto st = simData->ScalarData(entry);
    auto x = simData->X();
    // Set heights for particles in the window
    for (int i=0; i<simData->size_owned(); ++i) {
      if (x[i][0] < bound) st[i] = x[i][1]; // Set stripex to be the y value.
    }
  }

}