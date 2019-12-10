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
    RealType time = gflow->getElapsedTime();
    if (time-lastUpdate<updateDelay || temperature<=0) return;

    if (sim_dimensions>4) throw false;

    // Get data
    int size = simData->size();
    auto v = simData->V();
    auto f = simData->F();
    auto rd = simData->Sg();
    // Precomputed values, assumes Kb = 1
    RealType DT1 = temperature/(6.*viscosity*PI);
    RealType Df1 = sqrt(2.*DT1*(time-lastUpdate));
    // Add a random force
    RealType force[4], drag[4], strength;
    for (int n=0; n<size; ++n) {
      RealType Df2 = sqrt(1./rd(n));
      // Random normal direction
      randomNormalVec(force, sim_dimensions);
      // Random strength - 'temperature' is from the viscous medium
      strength = randNormal();
      scalarMultVec(Df1*Df2*strength, force, sim_dimensions);
      // Drag force - also from the viscous medium
      copyVec(v(n), drag, sim_dimensions);
      scalarMultVec(RealType(6.*PI*viscosity*rd(n)), drag, sim_dimensions);
      // Add total force
      plusEqVec (f(n), force, sim_dimensions);
      minusEqVec(f(n),  drag, sim_dimensions);
    }

    // Update time point - we shouldn't update this before now, so 
    // that (time-lastUpdate) will be correct (instead of 0)
    lastUpdate = time;
  }

}