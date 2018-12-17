#include "base.hpp"
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

  class GFlow* Base::getGFlow() {
    return gflow;
  }

  class SimData* Base::getSimData() {
    return simData;
  }

  class Integrator* Base::getIntegrator() {
    return integrator;
  }

  class DomainBase* Base::getDomain() {
    return domain;
  }

  class DataMaster* Base::getDataMaster() {
    return dataMaster;
  }

  class ForceMaster* Base::getForceMaster() {
    return forceMaster;
  }

}
