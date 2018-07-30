#include "domainbase.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : Base(gflow) {};

  const int* DomainBase::getDims() const {
    return dims;
  }

  const RealType* DomainBase::getWidths() const {
    return widths;
  }

  int DomainBase::getNumCells() const {
    int total = 1;
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) total *= dims[d];
    return total;
  }

  RealType DomainBase::getSkinDepth() const {
    return skinDepth;
  }

  RealType DomainBase::getCutoff() const {
    return cutoff;
  }

  int DomainBase::getNumberOfRemakes() const {
    return number_of_remakes;
  }

}