#include "death.hpp"
// Other files
#include "../base/integrator.hpp"

namespace GFlowSimulation {

  DeathRate::DeathRate(GFlow *gflow) : Modifier(gflow) {};

  void DeathRate::setRates(const vector<RealType>& rates) {
    deathRates = rates;
  }

  void DeathRate::pre_forces() {
    RealType dt = Base::integrator->getTimeStep();
    for (int i=0; i<Base::simData->number; ++i) {
      int type = Base::simData->type[i];
      if (deathRates.size()<=type) continue;
      RealType rate = deathRates[type];
      // Probability is 1-exp(-rate*dt) ~ rate*dt
      if (drand48()<rate*dt) Base::simData->markForRemoval(i);
    }
  }

}