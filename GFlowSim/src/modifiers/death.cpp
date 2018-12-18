#include "death.hpp"
// Other files
#include "../base/integrator.hpp"

namespace GFlowSimulation {

  DeathRate::DeathRate(GFlow *gflow) : Modifier(gflow) {};

  DeathRate::DeathRate(GFlow *gflow, const vector<RealType>& dr) : Modifier(gflow), deathRates(dr) {};

  void DeathRate::setRates(const vector<RealType>& rates) {
    deathRates = rates;
  }

  void DeathRate::pre_forces() {
    RealType dt = Base::integrator->getTimeStep();
    int size = Base::simData->size();
    for (int i=0; i<size; ++i) {
      int type = Base::simData->Type(i);
      if (deathRates.size()<=type || type<0) continue;
      RealType rate = deathRates[type];
      // Probability is 1-exp(-rate*dt) ~ rate*dt
      if (drand48()<rate*dt) Base::simData->markForRemoval(i);
    }
  }

}