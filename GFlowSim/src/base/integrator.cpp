#include "integrator.hpp"
// Other files
#include "simdata.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  Integrator::Integrator(GFlow *gflow) : Base(gflow), dt(0.0001), min_dt(1e-5), max_dt(0.001), target_steps(18), step_delay(0), step_count(step_delay) {};

  void Integrator::pre_integrate() {
    // Set step count so a check is triggered on the first step
    step_count = step_delay;
  }

  void Integrator::pre_step() {
    if (step_count < step_delay) {
      ++step_count;
      return;
    }
    // Call Base's pre_step (I don't think it does anything right now, but best to be safe)
    Base::pre_step();
    // Check the velocity components of all the particles
    RealType *v = simData->V_arr(), *sg = simData->Sg();
    // The minimum time a particle would take to traverse its own radius
    RealType minT = 1;

    // Find minT
    #if SIMD_TYPE==SIMD_NONE
    // Do serially
    for (int i=0; i<sim_dimensions*simData->number; ++i) {
      RealType mint = sg[i/sim_dimensions]/fabs(v[i]);
      if (mint<minT) minT = mint;
    }
    #else 
    // Do as much as we can in parallel
    simd_float MinT = simd_set1(1.);
    int i=0;
    for (; i<sim_dimensions*simData->number-simd_data_size; i += simd_data_size) {
      simd_float V = simd_abs(simd_load(&v[i]));
      simd_float Sg = simd_load_constant<DIMENSIONS>(sg, i);
      simd_float Mint = Sg / V;
      simd_float mask = simd_less_than(Mint, MinT);
      simd_update_masked(MinT, Mint, mask);
    }
    // Consolidate MinT
    for (int d=0; d<simd_data_size; ++d) {
      RealType mint = simd_get(d, MinT);
      if (mint<minT) minT = mint;
    }
    // Do the last part serially
    for (; i<sim_dimensions*simData->number; ++i) {
      RealType mint = sg[i/sim_dimensions]/fabs(v[i]);
      if (mint<minT) minT = mint;
    }
    #endif

    // Scale by a dimensional factor, since we just looked at every component of velocity separately
    minT /= sqrt(sim_dimensions);

    // Set the timestep
    dt = minT * 1./static_cast<RealType>(target_steps);
    if (dt>max_dt) dt = max_dt;
    else if (dt<min_dt) dt = min_dt;

    // Reset step count
    step_count = 0;
  }

  void Integrator::post_forces() {
    Base::post_forces();

    /*
    // Check the velocity components of all the particles
    RealType *f = simData->F_arr(), *im = simData->Im();

    RealType maxA = 0., aveA = 0;
    for (int i=0; i<sim_dimensions*simData->number; ++i) {
      RealType acc = im[i/sim_dimensions]*fabs(f[i]);
      if (acc>maxA) maxA = acc;

      aveA += acc;
    }

    RealType t = 1./maxA;
    if (t<dt) dt = t;
    if (dt<min_dt) dt = min_dt;

    cout << maxA << " " << aveA/(sim_dimensions*simData->number) << ", Suggest: " << t << endl;
    cout << dt << endl << endl;
    */
  }

  RealType Integrator::getTimeStep() {
    return dt;
  }

  void Integrator::setDT(RealType t) {
    dt = t;
  }

}