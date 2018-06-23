#include "force.hpp"

namespace GFlowSimulation {

  Force::Force(GFlow *gflow) : Base(gflow), typeMap(nullptr) {};

  Force::~Force() {
    if (typeMap)    delete [] typeMap;
  }

  void Force::clearVerletList() {
    verletList.clear();
  }

}