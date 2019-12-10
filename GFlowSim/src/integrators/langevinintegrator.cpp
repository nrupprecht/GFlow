#include "langevinintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../utility/printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LangevinIntegrator::LangevinIntegrator(GFlow *gflow) : LangevinTypeIntegrator(gflow, 0.025, DEFAULT_VISCOSITY) {};

  LangevinIntegrator::LangevinIntegrator(GFlow *gflow, RealType temp) : LangevinTypeIntegrator(gflow, temp, DEFAULT_VISCOSITY) {};

  void LangevinIntegrator::pre_forces() {
    // Start the timer
    timer.start();

    // Call parent class
    Integrator::pre_forces();

    // --- First half kick

    // Number of (real - non ghost) particles
    int total = simData->size_owned()*sim_dimensions;
    if (total==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    auto x = simData->X(), v = simData->V(), f = simData->F();
    auto im = simData->Im();

    // Update velocities
    for (int i=0; i<total; ++i) 
      v[i] += hdt*im[i/sim_dimensions]*f[i];

    // Update positions
    for (int i=0; i<total; ++i)
      x[i] += dt*v[i];

    // Stop timer
    timer.stop();
  }

  void LangevinIntegrator::post_forces() {
    // Start the timer
    timer.start();

    // Call to parent class
    Integrator::post_forces();
    
    // --- Second half kick

    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Number of (real - non ghost) particles
    int size = simData->size_owned();
    // Get arrays
    auto x = simData->X(), v = simData->V(), f = simData->F();
    auto im = simData->Im(), sg = simData->Sg();

    // Add random noise - we don't need to do this every time
    RealType time = Base::gflow->getElapsedTime();
    if (updateDelay<time-lastUpdate && temperature>0) {
      // Precomputed values, assumes Kb = 1
      RealType Df1 = sqrt(2.*drift1*(time-lastUpdate));
      // Add a random force to all spatial degrees of freedom
      for (int i=0; i<size*sim_dimensions; ++i) {
        RealType Df2 = sqrt(1./sg(i/sim_dimensions));
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
    }

    // Stop the timer
    timer.stop();
  }

}
