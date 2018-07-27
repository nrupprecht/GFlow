#include "hard_sphere.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Force(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {};

  
  void HardSphere::calculateForces() const {
    //verletList.forceLoop(this);
    //return;
    
    int nverlet = verletList.vlSize(), id1(0), id2(0); // List length, id pointers
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

    /*
    // --- Go through all particles
    for (int h=0; h<nheads-1; ++h) {
      h0 = heads[h]; 
      h1 = heads[h+1];    // This delimits the end of this part of the verlet list
      id1 = -verlet[h0++]; //--> // First particle head might interact with is the one after the head
      sigma = sg[id1];    // Interaction radius of the head particle
      for (; h0<h1; ++h0) {
        id2 = verlet[h0];
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
    // Last part of the lists - there is no "next head" to delimit the end, the end is the end of the list
    h0 = heads[nheads-1]; // Last head
    id1 = -verlet[h0++]; //--> // First particle is the one after the head
    sigma = sg[id1];
    for (; h0<nverlet; ++h0) {
      id2 = verlet[h0];
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sigma + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal);
        //normalVec(displacement, normal); // Get the normal vector
        // Calculate force strength
        forceStrength(F, normal, distance, id1, id2);
        // Add force
        plusEqVec (f[id1], F);
        minusEqVec(f[id2], F);        
      }
    }
    */
  }

/*
  void HardSphere::calculateForces() const {
    // Id pointers for the particles
    int id1(0), id2(0);
    // Set verlet list to begin
    if(!verletList.begin(id1)) return;
    // Get the data we need
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags
    // Get verlet list data
    RealType F[DIMENSIONS];
    // --- Go through all particles
    while (verletList.next(id1, id2)) {
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal); // Normalize distance -> normal
        // Calculate force strength
        forceStrength(F, normal, distance, id1, id2);
        // Add force
        plusEqVec (f[id1], F);
        minusEqVec(f[id2], F);
      }
    }
  }
  */

  void HardSphere::forceKernel(int id1, int id2) const {
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS]; // To calculate displacement, normal vector
    RealType F[DIMENSIONS];

    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // Get displacement
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

  void HardSphere::setRepulsion(RealType r) { 
    repulsion = r; 
  }

  inline void HardSphere::forceStrength(RealType *F, const RealType *normal, const RealType distance, const int id1, const int id2) const {
    scalarMultVec(repulsion*(simData->Sg(id1) + simData->Sg(id2) - distance), normal, F);
  }

}
