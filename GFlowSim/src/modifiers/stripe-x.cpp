#include "stripe-x.hpp"

namespace GFlowSimulation {

  StripeX::StripeX(GFlow *gflow) : Modifier(gflow), window(1.0), lastUpdate(0), updateDelay(0.1) {};

  void StripeX::pre_integrate() {
    // Get the strip x entry.
    entry = simData->requestScalarData("StripeX");
    // Entry must be positive. Otherwise, something is wrong.
    if (entry<0) throw false;
    // Initialize all values
    auto st = simData->ScalarData(entry);
    auto x = simData->X();
    // Set initial heights
    for (int i=0; i<simData->size_owned(); ++i)
      st(i) = x(i, 1);
  }
  
  void StripeX::post_forces() {
    // Get the time
    RealType time = Base::gflow->getElapsedTime();
    if (time-lastUpdate<updateDelay || entry<0) return;

    // Set particles in the window
    RealType bound = simData->getBounds().min[0] + window;
    auto x = simData->X();
    auto st = simData->ScalarData(entry);

    // Set heights for particles in the window
    for (int i=0; i<simData->size_owned(); ++i)
      if (x(i, 0)<bound) st(i) = x(i, 1);

    // Update time point
    lastUpdate = time;
  }

}