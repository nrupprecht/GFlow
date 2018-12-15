#include "integrator.hpp"
// Other files
#include "simdata.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  Integrator::Integrator(GFlow *gflow) : Base(gflow), dt(0.0001), adjust_dt(true), min_dt(1e-6), max_dt(0.002), target_steps(20), step_delay(10), step_count(step_delay+1) {};

  void Integrator::pre_integrate() {
    // Set step count so a check is triggered on the first step
    step_count = step_delay;

    // Compute average radius
    characteristic_length = 0;
    for (int n=0; n<simData->size(); ++n) {
      if (simData->Type(n)<0) continue;
      characteristic_length += simData->Sg(n);
    }
    characteristic_length /= static_cast<RealType>(simData->number());
  }

  void Integrator::pre_step() {
    // Call Base's pre_step (I don't think it does anything right now, but best to be safe)
    Base::pre_step();
    // If we are not adjusting dt, we are done.
    if (!adjust_dt) return;
    // Check if enough time has gone by
    if (step_count < step_delay) {
      ++step_count;
      return;
    }
    // Check the velocity components of all the particles
    RealType *v = simData->V_arr(), *sg = simData->Sg();
    const int total = sim_dimensions*simData->size();

    // Find maxV
    RealType maxV = 0;
    #if SIMD_TYPE==SIMD_NONE
    // Do serially
    for (int i=0; i<total; ++i)
      if (maxV<fabs(v[i])) maxV = fabs(v[i]);
    #else 
    // Do as much as we can in parallel
    simd_float MaxV = simd_set1(0.);
    int i=0;
    for (; i<total-simd_data_size; i += simd_data_size) {
      simd_float V = simd_abs(simd_load(&v[i]));
      simd_float mask = simd_less_than(MaxV, V);
      simd_update_masked(MaxV, V, mask);
    }
    // Consolidate MaxV
    for (int d=0; d<simd_data_size; ++d) {
      RealType mv = simd_get(d, MaxV);
      if (maxV<mv) maxV = mv;
    }
    // Do the last part serially
    for (; i<total; ++i)
      if (maxV<fabs(v[i])) maxV = fabs(v[i]);
    #endif

    // The minimum time any object takes to cover a characteristic length
    // @todo There should be a systematic finding of the number 0.05.
    RealType minT = characteristic_length/(maxV*sqrt(sim_dimensions));

    // Set the timestep
    dt = minT * 1./static_cast<RealType>(target_steps);
    if (dt>max_dt) dt = max_dt;
    else if (dt<min_dt) dt = min_dt;
    // Reset step count
    step_count = 0;
  }

  RealType Integrator::getTimeStep() {
    return dt;
  }

  void Integrator::setDT(RealType t) {
    dt = t;
  }

  void Integrator::setAdjustDT(bool d) {
    adjust_dt = d;
  }

  void Integrator::setTargetSteps(int s) {
    target_steps = max(1, s);
  }

  void Integrator::setStepDelay(int s) {
    step_delay = max(0, s);
  }

  void Integrator::setMaxDT(RealType t) {
    if (t>0) max_dt = t;
  }

  void Integrator::setMinDT(RealType t) {
    if (t>0) min_dt = t;
  }

}