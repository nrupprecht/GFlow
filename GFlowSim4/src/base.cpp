#include "base.hpp"
#include "gflow.hpp"
#include "simdata.hpp"
#include "sectorization.hpp"
#include "neighbors.hpp"
#include "communicator.hpp"

namespace GFlowSimulation {

  Base::Base(GFlow *gf) {
    gf->initializeBase(this);
  }

  void Base::initialize() {
    // Reset subpointers
    gflow->initializeBase(this);
  }

}