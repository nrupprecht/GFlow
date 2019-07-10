#include "domainbase.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : InteractionHandler(gflow) {
    // Allocate arrays
    dims     = new int[sim_dimensions];
    widths   = new RealType[sim_dimensions];
    inverseW = new RealType[sim_dimensions];
    // Set to zero
    zeroVec(dims, sim_dimensions);
    zeroVec(widths, sim_dimensions);
    zeroVec(inverseW, sim_dimensions);
  }; 

  DomainBase::~DomainBase() {
    if (dims)     delete [] dims;
    if (widths)   delete [] widths;
    if (inverseW) delete [] inverseW;
    dims     = nullptr;
    widths   = nullptr;
    inverseW = nullptr;
  }

  const int* DomainBase::getDims() const {
    return dims;
  }

  const RealType* DomainBase::getWidths() const {
    return widths;
  }

  int DomainBase::getNumCells() const {
    int total = 1;
    for (int d=0; d<sim_dimensions; ++d) total *= dims[d];
    return total;
  }

  RealType DomainBase::getCutoff() const {
    return min_small_cutoff;
  }

  void DomainBase::setSkinDepth(RealType s) {
    InteractionHandler::setSkinDepth(s);
    min_small_cutoff = 2*max_small_sigma + skin_depth;
  }

  void DomainBase::setBounds(const Bounds& bnds) {
    InteractionHandler::setBounds(bnds);
    // Recalculate dimensions. The domain must be initialized to be able to do this.
    if (initialized) calculate_domain_cell_dimensions();
  }

}
