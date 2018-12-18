#include "overdampedlangevin.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  OverdampedLangevinIntegrator::OverdampedLangevinIntegrator(GFlow *gflow) : Integrator(gflow), viscosity(DEFAULT_VISCOSITY), 
  temperature(10.), lastUpdate(0), updateDelay(0.01)
  {
    drift1 = temperature/(6.*viscosity*PI);
  }

  OverdampedLangevinIntegrator::OverdampedLangevinIntegrator(GFlow *gflow, RealType T) : Integrator(gflow), viscosity(DEFAULT_VISCOSITY), 
  temperature(T), lastUpdate(0), updateDelay(0.05)
  {
    drift1 = temperature/(6.*viscosity*PI);
  }

  void OverdampedLangevinIntegrator::post_forces() {
    // Call to parent class
    Integrator::post_forces();
    
    // Number of (real - non ghost) particles
    int size = simData->size();
    if (size==0) return;
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
        RealType strength = drand48()-0.5; // randNormal();
        f[i] += Df1*Df2*strength;
      }
      lastUpdate = time;
    }

    // Update positions (there are no velocities)
    for (int i=0; i<size*sim_dimensions; ++i) {
      int id = i/sim_dimensions;
      v[i] = viscosity*f[i]*im[id];
      x[i] += v[i]*Integrator::dt;
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(x[i]));
      assert(!isnan(v[i]));
      assert(fabs(v[i])<MAX_REASONABLE_V);
      #endif 
    }

  }

  void OverdampedLangevinIntegrator::setViscosity(RealType eta) {
    viscosity = eta;
    drift1 = temperature/(6.*viscosity*PI);
  }

  void OverdampedLangevinIntegrator::setTemperature(RealType T) {
    temperature = T;
    drift1 = temperature/(6.*viscosity*PI);
  }

}
