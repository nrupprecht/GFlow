#include "base.hpp"
#include "gflow.hpp"
#include "simdata.hpp"
#include "sectorization.hpp"
#include "neighbors.hpp"
#include "communicator.hpp"

namespace GFlowSimulation {

  Base::Base(GFlow *gflow) {
    this->gflow   = gflow;
    simData       = gflow->simData;
    sectorization = gflow->sectorization;
    neighbors     = gflow->neighbors;
    communicator  = gflow->communicator;
    integrator    = gflow->integrator;
  }

  void Base::initialize() {
    // Reset subpointers
    simData       = gflow->simData;
    sectorization = gflow->sectorization;
    neighbors     = gflow->neighbors;
    communicator  = gflow->communicator;
    integrator    = gflow->integrator;
  }

}