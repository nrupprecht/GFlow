#include "dataobject.hpp"

namespace GFlowSimulation {

  DataObject::DataObject(GFlow *gflow, string name) : Base(gflow), dataName(name) {};

  string DataObject::getName() {
    return dataName;
  }

}