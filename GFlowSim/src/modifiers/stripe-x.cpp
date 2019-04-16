#include "stripe-x.hpp"

namespace GFlowSimulation {

  StripeX::StripeX(GFlow *gflow) : Modifier(gflow), window(0.5), lastUpdate(0), updateDelay(0.1) {};

  void StripeX::pre_integrate() {
    // Get the strip x entry.
    entry = simData->requestScalarData("StripeX");
    // Entry must be positive. Otherwise, something is wrong.
    if (entry<0) throw false;
    // Initialize all values
    RealType *st = simData->ScalarData(entry), *x = simData->X_arr();
    // Set initial heights
    for (int i=0; i<simData->size()*sim_dimensions; i+=sim_dimensions)
      st[static_cast<int>(i/sim_dimensions)] = x[i+1];
  }
  
  void StripeX::post_forces() {
    // Get the time
    RealType time = Base::gflow->getElapsedTime();
    if (time-lastUpdate<updateDelay || entry<0) return;

    // Set particles in the window
    RealType bound = simData->getBounds().min[0] + window;
    RealType *st = simData->ScalarData(entry), *x = simData->X_arr();

    // Set heights for particles in the window
    for (int i=0; i<simData->size()*sim_dimensions; i+=sim_dimensions)
      if (x[i] < bound) st[static_cast<int>(i/sim_dimensions)] = x[i+1];

    // Update time point
    lastUpdate = time;
  }

}