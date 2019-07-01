#include "domainbase.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : InteractionHandler(gflow), domain_bounds(sim_dimensions), bounds(sim_dimensions) {
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

  void DomainBase::calculate_max_small_sigma() {
    // Make sure force master has interaction array set up
    forceMaster->initialize_does_interact();

    // Find average sigma
    RealType sigma = 0, max_sigma = 0;
    int count = 0;
    for (int n=0; n<Base::simData->size(); ++n) {
      // Check that the type is valid, and is an interacting type
      int type = Base::simData->Type(n);
      if (type<0 || !Base::forceMaster->typeInteracts(type)) 
        continue;
      // Get the cutoff radius, use in the calculation
      RealType s = Base::simData->Sg(n) * forceMaster->getMaxCutoff(type);
      sigma += s;
      if (s>max_sigma) max_sigma = s;
      ++count;
    }
    if (count>0) sigma /= count;
    else {
      sigma = Base::simData->Sg(0) * forceMaster->getMaxCutoff(simData->Type(0));
      max_sigma = sigma;
    }

    // Threshhold sigma is between average and maximum sigma
    RealType threshold = 0.5*(sigma + max_sigma), max_under = sigma;
    if (threshold!=sigma) {
      for (int n=0; n<Base::simData->size(); ++n) {
        // Check that the type is valid, and is an interacting type
        int type = Base::simData->Type(n);
        if (type<0 || !Base::forceMaster->typeInteracts(type)) 
          continue;
        // Get the cutoff radius, use in the calculation
        RealType s = Base::simData->Sg(n) * forceMaster->getMaxCutoff(simData->Type(n));
        if (s<threshold && max_under<s) max_under = s;
      }
    }
    max_small_sigma = 1.025*max_under;
  }

}
