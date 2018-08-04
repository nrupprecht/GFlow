#include "hardsphere-a.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  HardSphereA::HardSphereA(GFlow *gflow) : Force(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION), sigma(0.05) {};

  void HardSphereA::calculateForces() const {
    int nverlet = verletList.vlSize(), id1(0), id2(0); // List length, id pointers
    if (nverlet==0) return; // No forces to calculate

    // Get the data we need
    RealType **x = Base::simData->x, **f = Base::simData->f;
    int *type = Base::simData->type;
    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // Get verlet list data
    const int *verlet = verletList.getVerlet();
    RealType F[DIMENSIONS];

    // --- Go through all particles
    for (int i=0; i<nverlet; i+=2) {
      id1 = verlet[i];
      id2 = verlet[i+1];
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(2.*sigma)) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal);
        // Calculate force strength
        forceStrength(F, normal, distance, id1, id2);
      }
    }
  }

  //! @param r The desired repulsion.
  void HardSphereA::setRepulsion(RealType r) {
    repulsion = r;
  }

  //! @param s The desired cutoff radius.
  void HardSphereA::setSigma(RealType s) {
    sigma = s;
  }

  inline void HardSphereA::forceStrength(RealType *F, const RealType *normal, const RealType distance, const int id1, const int id2) const {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = Base::simData->type[id1]<0 ? 0 : 1.; //--
    RealType c2 = Base::simData->type[id2]<0 ? 0 : 1.; //--

    scalarMultVec(c1*c2*repulsion*(2.*sigma - distance), normal, F);
    // Add forces
    plusEqVec (Base::simData->f[id1], F);
    minusEqVec(Base::simData->f[id2], F);
  }

}