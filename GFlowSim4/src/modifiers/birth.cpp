#include "birth.hpp"
// Other files
#include "../base/integrator.hpp"

namespace GFlowSimulation {

  BirthRate::BirthRate(GFlow *gflow) : Modifier(gflow) {};

  void BirthRate::setRate(const vector<RealType>& rates) {
    birthRates = rates;
  }

  void BirthRate::pre_forces() {
    RealType dt = Base::integrator->getTimeStep();
    RealType X[DIMENSIONS], normal[DIMENSIONS];
    for (int i=0; i<Base::simData->number; ++i) {
      int type = Base::simData->type[i];
      if (birthRates.size()<=type) continue;
      RealType rate = birthRates[type];
      // Probability is 1-exp(-rate*dt) ~ rate*dt
      if (drand48()<rate*dt) {
        copyVec(Base::simData->x[i], X);
        randomNormalVec(normal);

        int global = Base::simData->getNextGlobalID();
        // Base::simData->addParticle();
      }
    }
  }

}