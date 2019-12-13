#include "overdampedintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/simd_utility.hpp"

template<int dimensions>
OverdampedIntegrator<dimensions>::OverdampedIntegrator(GFlow *gflow) : Integrator(gflow) {
  // Set dt
  dt = min_dt;
  // Set step delay to zero
  step_delay = 0;
};

template<int dimensions>
void OverdampedIntegrator<dimensions>::pre_integrate() {
  // Call this before
  Integrator::pre_integrate();
  // Set dt
  dt = min_dt;
}

template<int dimensions>
void OverdampedIntegrator<dimensions>::post_forces() {
  // Start the timer.
  timer.start();

  // Call to parent class.
  Integrator::post_forces();
  // Determine the propper time step.
  determine_time_step();

  // Number of (real - non ghost) particles.
  int total = dimensions*simData->size_owned();
  if (total==0) return;

  // Time step
  RealType dt = Integrator::dt;
  // Get arrays
  auto x = simData->X(), v = simData->V(), f = simData->F();
  auto im = simData->Im();

  // Update positions (there are no velocities).
  #if SIMD_TYPE==SIMD_NONE
  for (int i=0; i<total; ++i) {
    int id = i/sim_dimensions;
    x[i] += dampingConstant*im[id]*f[i]*dt;      
  }
  #else
  // Set simd vector with dampingConstant * dt.
  simd_float g_dt = simd_set1(dampingConstant*dt);
  int i;
  for (i=0; i<=total-simd_data_size; i+=simd_data_size) {
    simd_float _x = x.load_to_simd(i);
    simd_float _im = im.template valign_load_to_simd<dimensions>(i);
    simd_float _f = f.load_to_simd(i);
    // Do calculation.
    _x += g_dt * _im * _f;
    // Store result.
    x.store_simd(i, _x);
  }
  // Left overs.
  for (; i<total; ++i)
    x[i] += dampingConstant*im[i/dimensions]*f[i]*dt;
  #endif

  // Stop timer.
  timer.stop();
}

template<int dimensions>
void OverdampedIntegrator<dimensions>::setDamping(RealType d) {
  dampingConstant = d;
}

template<int dimensions>
inline void OverdampedIntegrator<dimensions>::determine_time_step() {
  // We need a condition on maximum force (actually acceleration), not maximum velocity, since dx/dt ~ F.
  if (!adjust_dt || simData->number()==0) return;

  // Check if enough time has gone by
  if (step_count < step_delay) {
    ++step_count;
    // If max acceleration is zero, it may because we have never checked what it is.
    maximum_acceleration = get_max_acceleration();
    // Return
    return;
  }
  // Reset step count
  step_count = 0;
  // Get the maximum acceleration of any particle
  maximum_acceleration = get_max_acceleration();
  // No data
  if (maximum_acceleration==0) return;
  // Set the timestep - acceleration is proportional to velocity (w/ constant dampingConstant)
  dt = characteristic_length/(dampingConstant*maximum_acceleration*static_cast<RealType>(target_steps)); 
  // Make sure timestep stays in the [min_dt, max_dt] bounds.
  if (dt>max_dt) dt = max_dt;
  else if (dt<min_dt) dt = min_dt;
}
