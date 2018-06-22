#include "dataobject.hpp"

namespace GFlowSimulation {

  DataObject::DataObject(GFlow *gflow, string name) : Base(gflow), dataName(name), delay(1./15.), lastRecording(-10.) {};

  string DataObject::getName() {
    return dataName;
  }

  void DataObject::setFPS(RealType fps) {
    delay = 1./fps;
  }

}
