#include "lennard_jones.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Force(gflow), strength(DEFAULT_LENNARD_JONES_STRENGTH), cutoff(DEFAULT_LENNARD_JONES_CUTOFF) {};

  void LennardJones::calculateForces() const {
    int nverlet = verletList.vlSize(), id1, id2; // List length, id pointers
    if (nverlet==0) return; // No forces to calculate

    // Get the data we need
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags
    RealType sigma; // Will hold the interaction radius of the head particle

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
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal);
        // Calculate force strength
        forceStrength(F, normal, distance, id1, id2);
        // Add force
        plusEqVec (f[id1], F);
        minusEqVec(f[id2], F);
      }
    }
  }

  void LennardJones::setStrength(RealType s) {
    strength = s;
  }

  inline void LennardJones::forceStrength(RealType *F, const RealType *normal, const RealType distance, const int id1, const int id2) const {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = Base::simData->type[id1]<0 ? 0 : 1.; //--
    RealType c2 = Base::simData->type[id2]<0 ? 0 : 1.; //--

    RealType *sg = simData->sg;
    // Calculate the magnitude
    RealType gamma = (sg[id1]+sg[id2]) / (distance*cutoff);
    RealType g3 = gamma*gamma*gamma, g6 = sqr(g3), g12 = sqr(g6);
    // The c1*c2 makes sure the magnitude is zero if either particles are of type -1
    RealType magnitude = c1*c2*24*strength*(2*g12 - g6)*(1./distance); 
    // Make vectorial
    scalarMultVec(magnitude, normal, F);
  }

}
