#include "integrator.hpp"
// Other files
#include "simdata.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  Integrator::Integrator(GFlow *gflow) : Base(gflow), dt(0.0001), adjust_dt(true), min_dt(1e-6), max_dt(0.002), 
    target_steps(20), step_delay(2), step_count(0), use_v(true), use_a(false) {};

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
    // Set dt to the minimum size
    if (adjust_dt) dt = min_dt;

    // cout << "Proposed max dt: " << sqrt(characteristic_length*characteristic_length/DEFAULT_HARD_SPHERE_REPULSION) << endl;
  }

  void Integrator::pre_step() {
    // If we are not adjusting dt, we are done.
    if (!adjust_dt) return;
    // Check if enough time has gone by
    if (step_count < step_delay) {
      ++step_count;
      return;
    }

    // Reset step count
    step_count = 0;
    // Get the maximum velocity
    RealType maxV = -1., maxA = -1., dt_v = 1., dt_a = 1.;
    if (use_v) {
      maxV = get_max_velocity();
      dt_v = characteristic_length*1./(maxV*static_cast<RealType>(target_steps));
    }
    if (use_a) {
      maxA = get_max_acceleration();
      dt_a = 10*sqrt(characteristic_length)*1./(maxA*static_cast<RealType>(target_steps));
    }
    // No information. Maybe this is the start of a run.
    if (maxV==0 && maxA==0) return;
    if (isnan(maxV) || isnan(maxA)) throw NanValue("Integrator pre-step detected NAN value.");
    // Set the timestep
    RealType dt_c = min(dt_v, dt_a); // Candidate dt
    dt = dt_c<dt ? dt_c : 0.9*dt + 0.1*dt_c;

    //cout << gflow->getElapsedTime() << ", V=" << maxV << ", A=" << maxA << "\t";

    if (dt>max_dt) dt = max_dt;
    else if (dt<min_dt) dt = min_dt;

    if (maxV>10) throw false;

    //cout << dt << " :: " << dt_v << ", " << dt_a << endl;
  }

  RealType Integrator::getTimeStep() {
    return dt;
  }

  void Integrator::setDT(RealType t) {
    dt = t;
  }

  void Integrator::setUseV(bool u) {
    use_v = u;
  }

  void Integrator::setUseA(bool u) {
    use_a = u;
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

  RealType Integrator::get_max_velocity() {
    // Check the velocity components of all the particles
    RealType *v = simData->V_arr();
    // Make sure the pointers are valid
    if (v==nullptr) return 0.;
    // Find maxV
    RealType maxV = 0;
    const int total = sim_dimensions*simData->size();

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

    // Return the max velocity
    return maxV*sqrt(sim_dimensions);
  }

  RealType Integrator::get_max_acceleration() {
    // Check the acceleration components of all the particles
    RealType *f = simData->F_arr(), *im = simData->Im();
    // Make sure the pointers are valid
    if (f==nullptr || im==nullptr) return 0.;
    // Reset the maximum acceleration of any particle.
    RealType maxA = 0.;
    const int total = sim_dimensions*simData->size();

    // Find minT
    #if SIMD_TYPE==SIMD_NONE
    // Do serially
    for (int i=0; i<total; ++i) {
      int id = i/sim_dimensions;
      RealType a = fabs(f[i]*im[id]);
      if (a>maxA) maxA = a;
    }
    #else 
    // Do as much as we can in parallel
    simd_float MaxA = simd_set1(0.);
    int i=0;
    for (; i<sim_dimensions*simData->size()-simd_data_size; i += simd_data_size) {
      simd_float F = simd_abs(simd_load(&f[i]));
      //simd_float Im = simd_load_constant(im, i, sim_dimensions);
      simd_float Im = simd_load_constant<2>(im, i);
      simd_float A = F*Im;
      simd_float mask = simd_less_than(MaxA, MaxA);
      simd_update_masked(MaxA, A, mask);
    }
    // Consolidate MaxA
    for (int d=0; d<simd_data_size; ++d) {
      RealType a = simd_get(d, MaxA);
      if (maxA<a) maxA = a;
    }
    // Do the last part serially
    for (; i<total; ++i) {
      RealType a = fabs(f[i]*im[i/sim_dimensions]);
      if (maxA<a) maxA = a;
    }
    #endif

    // Return the max acceleration
    return maxA*sqrt(sim_dimensions);
  }

}