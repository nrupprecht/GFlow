#include "TerminationCondition.hpp"

namespace GFlow {

  TerminationCondition::TerminationCondition(string m) : message(m) {};

  OutsideRegion::OutsideRegion() : TerminationCondition("A particle was outside the simulation bounds.") {};

  bool OutsideRegion::check(SimData* simData) {
    int domain_end = simData->getDomainEnd();
    Bounds simBounds = simData->getSimBounds();
    int *it = simData->getItPtr();
    for (int i=0; i<domain_end; ++i) {
      if (it[i]<0) continue;
      vec2 pos(simData->getPx(i), simData->getPy(i));
      if (!simBounds.contains(pos)) return true;
    }
    return false;
  }

  MaxVelocity::MaxVelocity(RealType v) : TerminationCondition("A particle's velocity exceeded the maximum allowed velocity."), maxVSqr(sqr(v)) {};

  bool MaxVelocity::check(SimData* simData) {
    int domain_end = simData->getDomainEnd();
    Bounds simBounds = simData->getSimBounds();
    int *it = simData->getItPtr();
    for (int i=0; i<domain_end; ++i) {
      if (it[i]<0) continue;
      RealType vsqr = sqr(simData->getVx(i)) + sqr(simData->getVy(i));
      if (vsqr>maxVSqr) return true;
    }
    return false;
  }

}
