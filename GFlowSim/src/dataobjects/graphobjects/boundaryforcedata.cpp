#include "boundaryforcedata.hpp"

namespace GFlowSimulation {

  BoundaryForceData::BoundaryForceData(GFlow *gflow) : GraphObject(gflow, "BDForce") {};

  void BoundaryForceData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // Store data. These functions work correctly with multiprocessor runs.
    gatherData(gflow->getElapsedTime(), gflow->getBoundaryForce());
  }

  RealType BoundaryForceData::getAverage() const {
    RealType force = 0;
    for (const auto bf : data)
      force += bf.second;
    return force/data.size();
  }

}