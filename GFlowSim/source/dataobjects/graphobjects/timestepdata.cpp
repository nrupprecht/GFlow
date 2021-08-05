#include <dataobjects/graphobjects/timestepdata.hpp>

using namespace GFlowSimulation;

TimeStepData::TimeStepData(GFlow *gflow)
    : GraphObject(gflow, "Timestep", "time", "time step") {};

void TimeStepData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  // Store data
  if (topology->getRank() == 0) {
    addEntry(Base::gflow->getElapsedTime(), gflow->getDT());
  }
}
