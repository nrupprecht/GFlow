#include "domaininteraction_test.hpp"
// Other files
#include "domain2d_test.hpp"

namespace GFlowSimulation {

  DomainInteraction::DomainInteraction(GFlow *gflow) : InteractionHandler(gflow) {};

  void DomainInteraction::executeKernel(Kernel, const RealType*, RealType*) const {
    Domain2d *domain2d = reinterpret_cast<Domain2d*>(Base::domain);

    // Make sure we are working with the right type of domain
    if (domain2d==nullptr) {
      cout << "This test class only works with domain2d. The domain is not of this type. Exiting.\n";
      exit(1);
    }

    // Have the domain calculate forces
    domain2d->calculateForces();
  }

}