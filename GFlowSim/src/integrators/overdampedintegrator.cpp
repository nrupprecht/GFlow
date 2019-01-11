#include "overdampedintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  OverdampedIntegrator::OverdampedIntegrator(GFlow *gflow) : Integrator(gflow), dampingConstant(DEFAULT_DAMPING_CONSTANT) {
    // Set dt
    dt = min_dt;
  };

  void OverdampedIntegrator::pre_integrate() {
    // Call this before
    Integrator::pre_integrate();
    // Set dt
    dt = min_dt;
  }

  void OverdampedIntegrator::pre_step() {
    // This overrides the default integrator prestep. The default integrator prestep adjusts the timestep so that
    // the fastest particle can only traverse its own radius in at least a set number of steps.
    // We need a condition on maximum force (actually acceleration), not maximum velocity, since dx/dt ~ F.
    if (!adjust_dt || simData->number()==0) return;
    // Check if enough time has gone by
    if (step_count < step_delay) {
      ++step_count;
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
    if (dt>max_dt) dt = max_dt;
    else if (dt<min_dt) dt = min_dt;
  }

  void OverdampedIntegrator::post_forces() {
    // Call to parent class
    Integrator::post_forces();

    // Number of (real - non ghost) particles
    int size = simData->size();
    if (size==0) return;

    // Time step
    RealType dt = Integrator::dt;
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im();

    // Update positions (there are no velocities)
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<size*sim_dimensions; ++i) {
      int id = i/sim_dimensions;
      x[i] += dampingConstant*im[id]*f[i]*dt;
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(f[i]));
      assert(!isnan(x[i]));
      assert(fabs(f[i])<MAX_REASONABLE_F);
      #endif 
    }
    #else
    // Get dampingConstant * dt
    simd_float g_dt = simd_set1(dampingConstant*dt);
    int i=0;
    for (; i<size*sim_dimensions-simd_data_size; i+=simd_data_size) {
      simd_float X = simd_load(&x[i]);
      simd_float _im;
      switch (sim_dimensions) {
        case 1: 
          _im = simd_load(&im[i]);
          break;
        case 2:
          _im = simd_load_constant<2>(im, i);
          break;
        case 3:
          _im = simd_load_constant<3>(im, i);
          break;
        case 4:
          _im = simd_load_constant<4>(im, i);
          break;
        default:
          throw false;
          break;
      }
      simd_float F = simd_load(&f[i]);
      simd_float dX = g_dt*_im*F;
      simd_float X_new  = X + dX; 
      simd_store(X_new, &x[i]);
    }
    // Left overs
    for (; i<size*sim_dimensions; ++i) {
      int id = i/sim_dimensions;
      x[i] += dampingConstant*im[id]*f[i]*dt;
    }
    #endif
  }

  void OverdampedIntegrator::setDamping(RealType d) {
    dampingConstant = d;
  }

}
