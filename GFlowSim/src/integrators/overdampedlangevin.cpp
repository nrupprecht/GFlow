#include "overdampedlangevin.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  OverdampedLangevinIntegrator::OverdampedLangevinIntegrator(GFlow *gflow) : LangevinTypeIntegrator(gflow, 1., DEFAULT_VISCOSITY) {};

  OverdampedLangevinIntegrator::OverdampedLangevinIntegrator(GFlow *gflow, RealType temp) : LangevinTypeIntegrator(gflow, temp, DEFAULT_VISCOSITY) {};

  void OverdampedLangevinIntegrator::post_forces() {
    // Start the timer
    timer.start();

    // Call to parent class
    Integrator::post_forces();
    
    // Number of (real - non ghost) particles
    int size = simData->size_owned();
    if (size==0) return;
    // Get arrays
    auto x = simData->X(), v = simData->V(), f = simData->F();
    auto im = simData->Im(), sg = simData->Sg();

    // Add random noise - we don't need to do this every time
    RealType time = Base::gflow->getElapsedTime(), dt = Integrator::dt;
    if (updateDelay<time-lastUpdate && temperature>0) {
      // Precomputed values, assumes Kb = 1
      RealType Df1 = sqrt(2.*drift1*(time-lastUpdate));
      // Add a random force to all spatial degrees of freedom
      for (int i=0; i<size*sim_dimensions; ++i) {
        RealType Df2 = sqrt(1./sg[i/sim_dimensions]);
        // Random strength - 'temperature' is from the viscous medium
        RealType strength = drand48()-0.5; // randNormal();
        f[i] += Df1*Df2*strength;
      }
      lastUpdate = time;
    }

    // Update positions (there are no velocities)
    for (int i=0; i<size*sim_dimensions; ++i) {
      int id = i/sim_dimensions;
      v[i] = viscosity*f[i]*im[id];
      x[i] += v[i]*dt;
    }

    // Stop timer
    timer.stop();
  }

}
