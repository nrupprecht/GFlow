#include "lennard_jones.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Force(gflow), strength(DEFAULT_LENNARD_JONES_STRENGTH), cutoff(DEFAULT_LENNARD_JONES_CUTOFF) {};

  void LennardJones::calculateForces() const {
    //verletList.forceLoop(this);
    //return;
    
    int nverlet = verletList.vlSize(), id1, id2; // List length, id pointers
    if (verletList.vlHSize()==0) return; // No forces to calculate

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
    for (int i=0; i<nverlet; ++i) {
      if (verlet[i]<0) {
        id1 = -verlet[i++]-1;
        sigma = sg[id1]; // Interaction radius of the head particle
      }
      id2 = verlet[i];
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sigma + sg[id2])) {
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

  void LennardJones::forceKernel(int id1, int id2) const {
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS]; // To calculate displacement, normal vector
    RealType F[DIMENSIONS];

    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
    // Check if the particles should interact
    RealType dsqr = sqr(displacement);
    if (dsqr < sqr(sg[id1] + sg[id2])) {
      RealType distance = sqrt(dsqr);
      scalarMultVec(1./distance, displacement);
      // Calculate force strength - displacement is now the normal vector
      forceStrength(F, displacement, distance, id1, id2);
      // Add force
      plusEqVec (f[id1], F);
      minusEqVec(f[id2], F);
    }
  }

  void LennardJones::setStrength(RealType s) {
    strength = s;
  }

  inline void LennardJones::forceStrength(RealType *F, const RealType *normal, const RealType distance, const int id1, const int id2) const {
    RealType *sg = simData->sg;
    // Calculate the magnitude
    RealType gamma = (sg[id1]+sg[id2]) / (distance*cutoff);
    RealType g3 = gamma*gamma*gamma, g6 = sqr(g3), g12 = sqr(g6);
    RealType magnitude = 24*strength*(2*g12 - g6)*(1./distance);
    // Make vectorial
    scalarMultVec(magnitude, normal, F);
  }

}
