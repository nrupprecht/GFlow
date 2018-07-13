#include "force.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  Force::Force(GFlow *gflow) : Base(gflow) {
    // Local copies of the data
    bounds = gflow->getBounds();
    copyVec(Base::gflow->getBCs(), boundaryConditions);
  };

  Force::~Force() {}

  void Force::initialize() {
    Base::initialize();
    // Local copies of the data
    bounds = gflow->getBounds();
    copyVec(Base::gflow->getBCs(), boundaryConditions);
  }

  int Force::lastHead() const {
    return verletList.lastHead();
  }

  int Force::vlSize() const {
    return verletList.vlSize();
  }

  int Force::vlHSize() const {
    return verletList.vlHSize();
  }

  const VerletList& Force::getVerletList() const {
    return verletList;
  }

  void Force::clearVerletList() {
    verletList.clear();
  }

  void Force::addVerletPair(int id1, int id2) {
    // Add the head if it is new
    if (id1!=lastHead()) verletList.addHead(id1);
    // Add id2 to the head
    verletList.addToHead(id2);
  }

  void Force::default_force_loop() {
    int nheads = verletList.vlHSize(), nverlet = verletList.vlSize();
    if (nheads==0) return; // No forces to calculate
    int h0, h1, id1, id2; // Head pointers, id pointers

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
          normalVec(displacement, normal); // Get the normal vector
          // Calculate force strength
          calculateStrength(F, normal, distance, id1, id2);
          // scalarMultVec(repulsion*(sigma + sg[id2] - distance), normal, F);
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
        // Calculate force strength
        calculateStrength(F, normal, distance, id1, id2);
        // scalarMultVec(repulsion*(sigma + sg[id2] - distance), normal, F);
        // Add force
        minusEqVec(f[id2], F);
        plusEqVec (f[id1], F);
      }
    }
  }

}