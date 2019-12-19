#include "integrator.hpp"
// Other files
#include "simdata.hpp"
#include "topology.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  Integrator::Integrator(GFlow *gflow) : Base(gflow), dt(0.0001), adjust_dt(true), min_dt(1e-6), max_dt(0.002), 
   target_steps(20), step_delay(10), step_count(0), use_v(true), use_a(false), characteristic_length(0.05) {};

  void Integrator::pre_integrate() {
    // Clear timer.
    clear_timer();
    
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
  }
  
  void Integrator::pre_step() {
    // If we are not adjusting dt, we are done.
    //! \todo Non-primary integrators should be allowed to request shorter timesteps.
    if (!adjust_dt || !isPrimaryIntegrator()) return;
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
      dt_v = characteristic_length/(maxV*static_cast<RealType>(target_steps));
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
    
    if (dt>max_dt) dt = max_dt;
    else if (dt<min_dt) dt = min_dt;

    #if USE_MPI == 1
      // Sync timesteps
      if (topology->getNumProc()>1) MPIObject::mpi_min(dt);
    #endif
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

  RealType Integrator::getMaxDT() const {
    return max_dt;
  }

  RealType Integrator::getMinDT() const {
    return min_dt;
  }

  RealType Integrator::get_max_velocity() const {
    // Check the velocity components of all the particles
    auto v = simData->V();

    // Find maxV
    RealType maxV = 0;

    #if SIMD_TYPE==SIMD_NONE
    // Do serially
    for (int i=0; i<simData->size_owned(); ++i) {
      real magnitude = magnitudeVec(v(i), sim_dimensions);
      if (maxV<magnitude) maxV = magnitude;
    }
    return maxV;
    #else 
    // Do as much as we can in parallel
    simd_float _maxv = simd_set1(0.);
    int i=0, total = simData->size_owned() * sim_dimensions;
    for (; i<total-simd_data_size; i += simd_data_size) {
      simd_float _va = simd_abs(v.load_to_simd(i)); // simd_abs(simd_load(&v[i]));
      simd_float _mask = simd_less_than(_maxv, _va);
      simd_update_masked(_maxv, _va, _mask);
    }
    // Consolidate _maxv
    for (int d=0; d<simd_data_size; ++d) {
      RealType mv = simd_get(d, _maxv);
      if (maxV<mv) maxV = mv;
    }
    // Do the last part serially
    for (; i<total; ++i)
      if (maxV<fabs(v[i])) maxV = fabs(v[i]);
    #endif
    // Return the max velocity
    return maxV*sqrt(sim_dimensions);
  }

  RealType Integrator::get_max_acceleration() const {
    // Check the acceleration components of all the particles
    auto f = simData->F();
    auto im = simData->Im();

    // Reset the maximum acceleration of any particle.
    RealType maxA = -1.;
    // Do serially
    for (int i=0; i<simData->size_owned(); ++i) {
      if (im[i]>0) {
        RealType a = magnitudeVec(f(i), sim_dimensions)*im[i];
        if (a>maxA) maxA = a;
      }
    }

    // Return the max acceleration
    return maxA;
  }

  bool Integrator::isPrimaryIntegrator() const {
    return (this==integrator);
  }

}
