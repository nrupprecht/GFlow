#include "dataobject.hpp"

namespace GFlowSimulation {

  DataObject::DataObject(GFlow *gflow, string name) : Base(gflow), dataName(name), delay(1./15.), lastRecording(-10.) {};

  string DataObject::getName() {
    return dataName;
  }

  void DataObject::setFPS(RealType fps) {
    delay = 1./fps;
  }

  bool DataObject::_check() {
    // Only record if enough time has gone by
    RealType time = Base::gflow->getElapsedTime();
    // If not enough time has gone by, return false
    if (time-lastRecording<delay) return false;
    // Otherwise, set last recording, return true
    lastRecording = time;
    return true;
  }

}
