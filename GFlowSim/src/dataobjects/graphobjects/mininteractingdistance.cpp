#include "mininteractingdistance.hpp"

namespace GFlowSimulation {

  MinInteractingDistance::MinInteractingDistance(GFlow* gflow) : GraphObject(gflow, "MinIntDist", "time", "min distance") {

    // This data object is currently unimplemented.
    throw false;

  };

  void MinInteractingDistance::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    
    /*
    RealType data_pack[] = { 10. }; // Random "large" number
    // Check each verlet list
    for (const auto it : gflow->getInteractions()) {
      // Pass in a null pointer for the param_pack, since we don't have any parameters
      it->executeKernel(&minimumDistanceKernel, nullptr, data_pack);
    }

    // Store data
    RealType time = Base::gflow->getElapsedTime();
    data.push_back(RPair(time, data_pack[0]));
    */
  }

}