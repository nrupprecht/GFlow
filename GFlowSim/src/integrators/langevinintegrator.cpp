#include "langevinintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../utility/printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LangevinIntegrator::LangevinIntegrator(GFlow *gflow) : Integrator(gflow), viscosity(DEFAULT_VISCOSITY), 
    temperature(0.025), lastUpdate(0), updateDelay(DEFAULT_TEMPERATURE_UPDATE_DELAY)
  {
    drift1 = temperature/(6.*viscosity*PI);
  }

  LangevinIntegrator::LangevinIntegrator(GFlow *gflow, RealType T) : Integrator(gflow), viscosity(DEFAULT_VISCOSITY), 
  temperature(T), lastUpdate(0), updateDelay(DEFAULT_TEMPERATURE_UPDATE_DELAY)
  {
    drift1 = temperature/(6.*viscosity*PI);
  }

  void LangevinIntegrator::pre_forces() {
    // --- First half kick

    // Number of (real - non ghost) particles
    int number = simData->number;
    if (number==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im;

    // Update velocities
    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
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
    for (int i=0; i<number*DIMENSIONS; ++i)
      x[i] += dt*v[i];
  }

  void LangevinIntegrator::post_forces() {
    // --- Second half kick

    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Number of (real - non ghost) particles
    int number = simData->number;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im, *sg = simData->sg;

    // Add random noise - we don't need to do this every time
    RealType time = Base::gflow->getElapsedTime();
    if (updateDelay<time-lastUpdate && temperature>0) {
      // Precomputed values, assumes Kb = 1
      RealType Df1 = sqrt(2.*drift1*(time-lastUpdate));
      // Add a random force to all spatial degrees of freedom
      for (int i=0; i<number*DIMENSIONS; ++i) {
        int id = i/DIMENSIONS;
        RealType Df2 = sqrt(1./sg[id]);
        // Random strength - 'temperature' is from the viscous medium
        RealType strength = drand48()-0.5; //randNormal();
        f[i] += Df1*Df2*strength;
      }
      lastUpdate = time;
    }

    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
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

  void LangevinIntegrator::setViscosity(RealType eta) {
    viscosity = eta;
    drift1 = temperature/(6.*viscosity*PI);
  }

  void LangevinIntegrator::setTemperature(RealType T) {
    temperature = T;
    drift1 = temperature/(6.*viscosity*PI);
  }

}
