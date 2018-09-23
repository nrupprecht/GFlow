#include "harmonicbond.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../base/simdata.hpp"
#include "../utility/printingutility.hpp" // For debugging

namespace GFlowSimulation {

  HarmonicBond::HarmonicBond(GFlow *gflow, int i1, int i2) : Modifier(gflow), id1(i1), id2(i2), springK(DEFAULT_SPRING_K) {
    // Use their current distance as the equilibrium distance
    RealType displacement[DIMENSIONS];
    subtractVec(Base::simData->x[i1], Base::simData->x[i2], displacement);
    eqDist = sqrt(sqr(displacement));
  };

  void HarmonicBond::post_forces() {
    // Get positions
    RealType *x1 = Base::simData->x[id1], *x2 = Base::simData->x[id2];
    // Figure out force
    RealType force[DIMENSIONS];
    getDisplacement(x1, x2, force, Base::gflow->getBounds(), Base::gflow->getBCs()); // Displacement
    // Distance
    RealType distance = sqrt(sqr(force));
    // Find x-x_eq
    RealType diffX = distance - eqDist;
    // Normalize [force]
    normalizeVec(force);
    // Find force
    scalarMultVec(springK*diffX, force); // F = - K \vec{x-x_eq}
    // Add the forces
    minusEqVec(Base::simData->f[id1], force);
    plusEqVec(Base::simData->f[id2], force);
  }

  void HarmonicBond::setSpringK(RealType sk) {
    springK = sk;
  }

}
