#include "overdampedintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/simd_utility.hpp"

template<int dimensions>
OverdampedIntegrator<dimensions>::OverdampedIntegrator(GFlow *gflow) : OverdampedIntegratorBase(gflow) {
  // Set dt
  dt = min_dt;
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

  // Potentially recalculate time step.
  if (adjust_dt && step_delay<=step_count) {
    calculate_time_step();
    step_count = 0;
  }
  else {
    ++step_count;
    // Recheck the maximum acceleration if it is zero (zero likely means that no data has been gathered).
    if (maximum_acceleration==0) calculate_time_step();
  }

  // Number of (real - non ghost) particles.
  int total = dimensions*simData->size_owned();
  if (total==0) return;

  // Time step.
  RealType dt = Integrator::dt;
  // Get arrays.
  auto x = simData->X(), v = simData->V(), f = simData->F();
  auto im = simData->Im();

  // Update positions (there are no velocities).
  #if SIMD_TYPE==SIMD_NONE
  for (int i=0; i<total; ++i) {
    int id = i/sim_dimensions;
    x[i] += dampingConstant*im[id]*f[i]*dt;      
  }
  #else
  // Set simd vector with (dampingConstant * dt).
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
  // Left overs
  for (; i<total; ++i)
    x[i] += dampingConstant*im[i/dimensions]*f[i]*dt;
  #endif

  // Stop timer.
  timer.stop();
}

template<int dimensions>
void OverdampedIntegrator<dimensions>::calculate_time_step() {
  // Get the maximum acceleration of any particle. 
  // We need a condition on maximum force (actually acceleration), not maximum velocity, since dx/dt ~ F.
  maximum_acceleration = get_max_acceleration();
  // If no data, just return.
  if (maximum_acceleration==0) return;
  // Set the timestep - acceleration is proportional to velocity (w/ constant dampingConstant)
  dt = characteristic_length/(dampingConstant*maximum_acceleration*static_cast<RealType>(target_steps));
  if (dt>max_dt) dt = max_dt;
  else if (dt<min_dt) dt = min_dt;

  #if USE_MPI == 1
    // Sync timesteps
    if (topology->getNumProc()>1) MPIObject::mpi_min(dt);
  #endif
}
