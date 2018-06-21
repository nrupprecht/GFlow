#include "force.hpp"

namespace GFlowSimulation {

  Force::Force(GFlow *gflow) : Base(gflow), verletList(nullptr), typeMap(nullptr) {};

  Force::~Force() {
    if (verletList!=nullptr) delete [] verletList;
    if (typeMap!=nullptr)    delete [] typeMap;
  }

}