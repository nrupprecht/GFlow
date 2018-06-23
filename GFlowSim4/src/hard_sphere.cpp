#include "hard_sphere.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Force(gflow), strength(DEFAULT_HARD_SPHERE_STRENGTH) {};

  void HardSphere::initialize() {
    // Set up verlet list
    
  }

  void HardSphere::calculateForces() {
    int nheads = verletList.vlHSize(), nverlet = verletList.vlSize();
    int h0, h1, id1, id2; // Head pointers, id pointers

    // Get the data we need
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS]; // To calculate displacement
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    bool wrap[DIMENSIONS]; 
    copyVec(Base::gflow->getWrap(), wrap); // Keep a local copy of the wrap frags
    RealType sigma; // Will hold the interaction radius of the head particle

    // Get verlet list data
    const int *verlet = verletList.getVerlet();
    const int *heads  = verletList.getHeads();
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
        getDisplacement(x[id1], x[id2], displacement, bounds, wrap);
        // Check if the particles should interact
        RealType dsqr = sqr(displacement);

        // cout << id1 << ", " << id2 << ": " << sqrt(dsqr) << " " << sigma + sg[id2] << endl;

        if (dsqr < sqr(sigma + sg[id2])) {
          RealType distance = sqrt(dsqr);
          scalarMultVec(strength*(sigma + sg[id2] - distance), displacement, F);
          // Add force
          minusEqVec(f[id1], F);
          plusEqVec (f[id2], F);
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
      getDisplacement(x[id1], x[id2], displacement, bounds, wrap);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sigma + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(100*(sigma + sg[id2] - distance), displacement, F);
        // Add force
        minusEqVec(f[id1], F);
        plusEqVec (f[id2], F);
      }
    }
  }

}