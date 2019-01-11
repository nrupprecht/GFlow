#include "langevinintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../utility/printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LangevinIntegrator::LangevinIntegrator(GFlow *gflow) : LangevinTypeIntegrator(gflow, 0.025, DEFAULT_VISCOSITY) {};

  LangevinIntegrator::LangevinIntegrator(GFlow *gflow, RealType temp) : LangevinTypeIntegrator(gflow, temp, DEFAULT_VISCOSITY) {};

  void LangevinIntegrator::pre_forces() {
    // --- First half kick

    // Number of (real - non ghost) particles
    int size = simData->size();
    if (size==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im();

    // Update velocities
    for (int i=0; i<size*sim_dimensions; ++i) {
      int id = i/sim_dimensions;
      v[i] += hdt*im[id]*f[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(f[i]));
      assert(!isnan(v[i]));
      assert(fabs(v[i])<MAX_REASONABLE_V);
      assert(fabs(f[i])<MAX_REASONABLE_F);
      #endif 
    }

    // Update positions
    for (int i=0; i<size*sim_dimensions; ++i)
      x[i] += dt*v[i];
  }

  void LangevinIntegrator::post_forces() {
    // Call to parent class
    Integrator::post_forces();
    
    // --- Second half kick

    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Number of (real - non ghost) particles
    int size = simData->size();
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im(), *sg = simData->Sg();

    // Add random noise - we don't need to do this every time
    RealType time = Base::gflow->getElapsedTime();
    if (updateDelay<time-lastUpdate && temperature>0) {
      // Precomputed values, assumes Kb = 1
      RealType Df1 = sqrt(2.*drift1*(time-lastUpdate));
      // Add a random force to all spatial degrees of freedom
      for (int i=0; i<size*sim_dimensions; ++i) {
        int id = i/sim_dimensions;
        RealType Df2 = sqrt(1./sg[id]);
        // Random strength - 'temperature' is from the viscous medium
        RealType strength = drand48()-0.5; //randNormal();
        f[i] += Df1*Df2*strength;
      }
      lastUpdate = time;
    }

    for (int i=0; i<size*sim_dimensions; ++i) {
      int id = i/sim_dimensions;
      // Drag force
      f[i] -= 6.*PI*viscosity*sg[id]*v[i];
      // Update velocity
      v[i] += hdt*im[id]*f[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(f[i]));
      assert(!isnan(v[i]));
      assert(fabs(v[i])<MAX_REASONABLE_V);
      assert(fabs(f[i])<MAX_REASONABLE_F);
      #endif 
    }
  }

}
