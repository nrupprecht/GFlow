#include "base.hpp"
// Other files.
#include "../gflow.hpp"

namespace GFlowSimulation {

  Base::Base(GFlow *gf) {
    if (gf) gf->initializeBase(this);
  }

  void Base::initialize() {
    // Reset subpointers
    if (gflow) gflow->initializeBase(this);
    else throw UnexpectedNullPointer("GFlow should not be null.");
  }

  class GFlow* Base::getGFlow() const {
    return gflow;
  }

  class SimData* Base::getSimData() const {
    return simData;
  }

  class Integrator* Base::getIntegrator() const {
    return integrator;
  }

  class InteractionHandler* Base::getHandler() const {
    return handler;
  }

  class DataMaster* Base::getDataMaster() const {
    return dataMaster;
  }

  class ForceMaster* Base::getForceMaster() const {
    return forceMaster;
  }

  class Topology* Base::getTopology() const {
    return topology;
  }

  int Base::getSimDimensions() const {
    return sim_dimensions;
  }

}
