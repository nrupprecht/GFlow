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

  int Base::memory_usage() {
    return 10*sizeof(int*);
  }

  GFlow* Base::getGFlow() const {
    return gflow;
  }

  shared_ptr<SimData> Base::getSimData() const {
    return simData;
  }

  Integrator* Base::getIntegrator() const {
    return integrator;
  }

  InteractionHandler* Base::getHandler() const {
    return handler;
  }

  DataMaster* Base::getDataMaster() const {
    return dataMaster;
  }

  ForceMaster* Base::getForceMaster() const {
    return forceMaster;
  }

  Topology* Base::getTopology() const {
    return topology;
  }

  int Base::getSimDimensions() const {
    return sim_dimensions;
  }

}
