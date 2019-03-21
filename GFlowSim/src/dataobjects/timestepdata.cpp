#include "timestepdata.hpp"

namespace GFlowSimulation {

  TimeStepData::TimeStepData(GFlow *gflow) : GraphObject(gflow, "Timestep", "time", "time step") {};

  void TimeStepData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Store data
    RealType time = Base::gflow->getElapsedTime();
    data.push_back(RPair(time, gflow->getDT()));
  }

}