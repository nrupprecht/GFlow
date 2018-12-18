#include "stripe-x.hpp"

namespace GFlowSimulation {

  StripeX::StripeX(GFlow *gflow) : Modifier(gflow), window(0.5), lastUpdate(0), updateDelay(0.1) {};

  void StripeX::pre_integrate() {
    entry = simData->request_scalar_data("StripeX");
    // Initialize all values
    RealType *st = simData->ScalarData(entry);
    for (int n=0; n<simData->size(); ++n)
      st[n] = simData->X(n, 1);
  }

  void StripeX::post_forces() {
    // Get the time
    RealType time = Base::gflow->getElapsedTime();
    if (time-lastUpdate<updateDelay || entry<0) return;
    // Set particles in the window
    RealType bound = simData->getBounds().min[0] + window;
    RealType *st = simData->ScalarData(entry);
    for (int n=0; n<simData->size(); ++n) {
      if( simData->X(n, 0) < bound) st[n] = simData->X(n, 1);
    }

    // Update time point
    lastUpdate = time;
  }

}