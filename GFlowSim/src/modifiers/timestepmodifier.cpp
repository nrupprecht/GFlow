#include "timestepmodifier.hpp"

namespace GFlowSimulation {

  TimestepModifier::TimestepModifier(GFlow *gflow) : Modifier(gflow), minSigma(1.0), vtollerance(0.07), atollerance(1), 
    lastCheck(0), delay(0.025), max_dt(DEFAULT_MAX_DT), min_dt(DEFAULT_MIN_DT) {}

  void TimestepModifier::pre_integrate() {
    SimData *simData = Base::simData;
    // Check
    if (simData==nullptr || simData->number==0) return;
    // Look for minimum sigma
    RealType *sg = simData->Sg();
    minSigma = sg[0];
    for (int n=0; n<simData->number; ++n)
      if (sg[n]<minSigma) minSigma = sg[n];
  }

  void TimestepModifier::post_forces() {
    // Make sure enough time has gone by
    RealType time = gflow->getElapsedTime();
    if (time-lastCheck<delay) return;
    lastCheck = time;

    // --- Find max quantities
    RealType maxV(0), maxA(0);
    SimData *simData = Base::simData;
    // Find maximum velocity
    for (int n=0; n<simData->number; ++n) {
      RealType v = magnitudeVec(simData->V(n));
      if (v>maxV) maxV = v;
    }
    // Find maximum acceleration
    RealType *im = simData->Im();
    for (int n=0; n<simData->number; ++n) {
      RealType a = magnitudeVec(simData->F(n))*im[n];
      if (a>maxA) maxA = a;
    }

    // --> Should really look for maximum dF/dx * v for particles

    // Propose timesteps
    RealType dt_v = vtollerance*minSigma/maxV, dt_a = atollerance/maxA;
    // Set the timestep
    RealType proposed_dt = min(dt_v, dt_a), current_dt = Base::gflow->getDT();
    proposed_dt = min(proposed_dt, max_dt);
    proposed_dt = max(min_dt, proposed_dt);
    Base::gflow->setDT(delay*(proposed_dt - current_dt) + current_dt);
  }

}