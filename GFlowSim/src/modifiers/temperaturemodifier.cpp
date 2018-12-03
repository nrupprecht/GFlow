#include "temperaturemodifier.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  TemperatureModifier::TemperatureModifier(GFlow *gflow, RealType t) : Modifier(gflow), 
    temperature(t), viscosity(DEFAULT_VISCOSITY), lastUpdate(0), updateDelay(0.01) {};

  TemperatureModifier::TemperatureModifier(GFlow *gflow, RealType t, RealType h) : Modifier(gflow), 
    temperature(t), viscosity(h), lastUpdate(0), updateDelay(0.01) {};

  void TemperatureModifier::post_forces() {
    // Get the time
    RealType time = Base::gflow->getElapsedTime();
    if (time-lastUpdate>updateDelay || temperature<=0) return;
    // Get data
    int number = Base::simData->number;
    RealType **v = Base::simData->V();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();
    // Precomputed values, assumes Kb = 1
    RealType DT1 = temperature/(6.*viscosity*PI);
    RealType Df1 = sqrt(2.*DT1*(time-lastUpdate));
    // Add a random force
    RealType *force = new RealType[sim_dimensions], *drag = new RealType[sim_dimensions], strength;
    for (int n=0; n<number; ++n) {
      RealType Df2 = sqrt(1./sg[n]);
      // Random normal direction
      randomNormalVec(force, sim_dimensions);
      // Random strength - 'temperature' is from the viscous medium
      strength = randNormal();
      scalarMultVec(Df1*Df2*strength, force, sim_dimensions);
      // Drag force - also from the viscous medium
      copyVec(v[n], drag, sim_dimensions);
      scalarMultVec(6.*PI*viscosity*sg[n], drag, sim_dimensions);
      // Add total force
      plusEqVec (f[n], force, sim_dimensions);
      minusEqVec(f[n],  drag, sim_dimensions);
    }

    // Update time point - we shouldn't update this before now, so 
    // that (time-lastUpdate) will be correct (instead of 0)
    lastUpdate = time;
    // Clean up
    delete [] force;
    delete [] drag;
  }

}