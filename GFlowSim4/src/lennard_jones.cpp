#include "lennard_jones.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Force(gflow), strength(DEFAULT_LENNARD_JONES_STRENGTH), cutoff(DEFAULT_LENNARD_JONES_CUTOFF) {};

  void LennardJones::calculateForces() {
    int nheads = verletList.vlHSize(), nverlet = verletList.vlSize();
    if (nheads==0) return; // No forces to calculate
    int h0, h1, id1, id2; // Head pointers, id pointers

    // Get the data we need
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    RealType sigma; // Will hold the interaction radius of the head particle

    // Get verlet list data
    const int *verlet = verletList.getVerlet();
    const int *heads  = verletList.getHeads();
    // Hold the forces in the same array
    RealType F[DIMENSIONS];
    // --- Go through all particles
    for (int h=0; h<nheads-1; ++h) {
      h0 = heads[h]; 
      h1 = heads[h+1];    // This delimits the end of this part of the verlet list
      id1 = verlet[h0++]; // First particle head might interact with is the one after the head
      sigma = sg[id1];    // Interaction radius of the head particle
      for ( ; h0<h1; ++h0) {
        id2 = verlet[h0];
        // Get the displacement between the particles
        getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
        // Check if the particles should interact
        RealType dsqr = sqr(displacement);
        if (dsqr < sqr(sigma + sg[id2])) {
          RealType distance = sqrt(dsqr);
          // Get the normal vector
          scalarMultVec(1./distance, displacement, normal);
          // normalVec(displacement, normal); 
          // Calculate force strength
          forceStrength(F, normal, distance, id1, id2);
          // Add force
          plusEqVec (f[id1], F);
          minusEqVec(f[id2], F);
        }
      }
    }
    // Last part of the lists - there is no "next head" to delimit the end, the end is the end of the list
    h0 = heads[nheads-1]; // Last head
    id1 = verlet[h0++];   // First particle is the one after the head
    sigma = sg[id1];
    for (; h0<nverlet; ++h0) {
      id2 = verlet[h0];
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sigma + sg[id2])) {
        RealType distance = sqrt(dsqr);
        normalVec(displacement, normal); // Get the normal vector
        // Calculate force strength
        forceStrength(F, normal, distance, id1, id2);
        // Add force
        plusEqVec (f[id1], F);
        minusEqVec(f[id2], F);
      }
    }
  }

  void LennardJones::forceKernel(int id1, int id2) {
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

  inline void LennardJones::forceStrength(RealType* F, RealType *normal, RealType distance, int id1, int id2) {
    RealType *sg = simData->sg;
    // Calculate the magnitude
    RealType gamma = (sg[id1]+sg[id2]) / (distance*cutoff);
    RealType g3 = gamma*gamma*gamma, g6 = sqr(g3), g12 = sqr(g6);
    RealType magnitude = 24*strength*(2*g12 - g6)*(1./distance);
    // Make vectorial
    scalarMultVec(magnitude, normal, F);
  }

}
