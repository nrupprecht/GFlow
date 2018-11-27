#include "overdampedintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  OverdampedIntegrator::OverdampedIntegrator(GFlow *gflow) : Integrator(gflow), dampingConstant(10*DEFAULT_DAMPING_CONSTANT) {
    max_dt = 0.001;
  };

  void OverdampedIntegrator::pre_step() {
    // Call Base's pre_step
    Base::pre_step();

    if (!adjust_dt) return;
    // Check if enough time has gone by
    if (step_count < step_delay) {
      ++step_count;
      return;
    }

    // This overrides the default integrator prestep. The default integrator prestep adjusts the timestep so that
    // the fastest particle can only traverse its own radius in at least a set number of steps.
    // We need a condition on maximum force (actually acceleration), not maximum velocity, since dx/dt ~ F.
    
    // Check the acceleration components of all the particles
    RealType *f = simData->F_arr(), *im = simData->Im(), *sg = simData->Sg();
    // The minimum time a particle would take to traverse its own radius
    //!  @todo A more nuanced thing to check would be how long it takes the fastest
    //!  particle to traverse the smallest radius, or the smallest radius "near" it.
    //!  The smallest radius in each subdivision could be found by binning.
    RealType minT = 1.; // Starting value

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
      simd_float F = simd_abs(simd_load(&f[i]));
      simd_float Im = simd_load_constant<DIMENSIONS>(im, i);
      simd_float Sg = simd_load_constant<DIMENSIONS>(sg, i);
      simd_float Mint = Sg / (F * Im);
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
      RealType mint = sg[i/sim_dimensions] / fabs(f[i]*im[i]);
      if (mint<minT) minT = mint;
    }
    #endif

    // Scale by a dimensional factor, since we just looked at every component of velocity separately
    minT /= sqrt(sim_dimensions);

    // Set the timestep
    dt = minT * 1./(static_cast<RealType>(target_steps) * dampingConstant);
    if (dt>max_dt) dt = max_dt;
    else if (dt<min_dt) dt = min_dt;

    // Reset step count
    step_count = 0;
  }

  void OverdampedIntegrator::post_forces() {
    // Call to parent class
    Integrator::post_forces();
    
    // Number of (real - non ghost) particles
    int number = simData->number;
    if (number==0) return;

    // Time step
    RealType dt = Integrator::dt;
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im();

    // Update positions (there are no velocities)
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
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
    for (; i<number*DIMENSIONS-simd_data_size; i+=simd_data_size) {
      simd_float X = simd_load(&x[i]);
      simd_float _im = simd_load_constant<DIMENSIONS>(im, i);
      simd_float F = simd_load(&f[i]);
      simd_float dX = g_dt*_im*F;
      simd_float X_new  = X + dX; 
      simd_store(X_new, &x[i]);
    }
    // Left overs
    for (; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      x[i] += dampingConstant*im[id]*f[i]*dt;
    }
    #endif
  }

  void OverdampedIntegrator::setDamping(RealType d) {
    dampingConstant = d;
  }

}