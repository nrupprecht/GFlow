#include "harmonicbond.hpp"
// Other files
#include "vectormath.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  HarmonicBond::HarmonicBond(GFlow *gflow, int i1, int i2) : Modifier(gflow), id1(i1), id2(i2), springK(DEFAULT_SPRING_K) {};

  void HarmonicBond::post_forces() {
    // Get positions
    RealType *x1 = Base::simData->x[id1], *x2 = Base::simData->x[id2];
    // Figure out force
    RealType force[DIMENSIONS];
    subtractVec(x1, x2, force); // Displacement
    scalarMultVec(springK, force); // F = - K \vec{x}
    // Add the force -- *** NOT SURE IF +/- SHOULD BE SWITCHED
    plusEqVec(Base::simData->f[id1], force);
    minusEqVec(Base::simData->f[id1], force);
  }

}
