#include "velocitylimiter.hpp"

namespace GFlowSimulation {

  VelocityLimiter::VelocityLimiter(GFlow *gflow, RealType maxv) : Modifier(gflow), maxV(maxv), lastUpdate(-1), updateDelay(0.01) {};

  void VelocityLimiter::post_forces() {
    // Get the time
    RealType time = Base::gflow->getElapsedTime();
    if (time-lastUpdate<updateDelay) return;
    // Get data
    int number = Base::simData->number;
    RealType **v = Base::simData->V();
    RealType maxVsqr = sqr(maxV);

    RealType *v_arr = Base::simData->V_arr();
    for (int i=0; i<number*sim_dimensions; ++i) {
      if (v_arr[i]>maxV) {
        Base::simData->markForRemoval(i);
        i = sim_dimensions*(i/sim_dimensions + 1);
      }
    }

    // Update record
    lastUpdate = time;
  }

}