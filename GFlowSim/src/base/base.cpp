#include "base.hpp"
#include "../gflow.hpp"

namespace GFlowSimulation {

  Base::Base(GFlow *gf) : sim_dimensions(gf->sim_dimensions) {
    if (gf) gf->initializeBase(this);
  }

  void Base::initialize() {
    // Reset subpointers
    gflow->initializeBase(this);
  }

}