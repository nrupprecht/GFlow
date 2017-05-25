#include "TerminationCondition.hpp"

namespace GFlow {

  TerminationCondition::TerminationCondition(string m) : message(m) {};

  OutsideRegion::OutsideRegion() : TerminationCondition("A particle was outside the simulation bounds") {};

  bool OutsideRegion::check(SimData* simData) {
    int domain_size = simData->getDomainSize();
    Bounds simBounds = simData->getSimBounds();
    for (int i=0; i<domain_size; ++i) {
      vec2 pos(simData->getPx(i), simData->getPy(i));
      if (!simBounds.contains(pos)) return true;
    }
    return false;
  }

}
