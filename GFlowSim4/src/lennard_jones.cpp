#include "lennard_jones.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Force(gflow), strength(DEFAULT_LENNARD_JONES_STRENGTH), cutoff(DEFAULT_LENNARD_JONES_CUTOFF) {};

  void LennardJones::calculateForces() const {
    // Check if there are forces to calculate. Virial is reset.
    if (!Force::initCalculation()) return; 

    // Get the data we need
    int nverlet = verletList.vlSize(), id1, id2; // List length, id pointers
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
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
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal);
        // Calculate force strength
        forceStrength(normal, distance, id1, id2);
      }
    }
  }

  void LennardJones::setStrength(RealType s) {
    strength = s;
  }

  //! \param[in] normal The normal of the displacement between the particles. Points from id2 to id1 (?). It is safe
  //! to change this parameter.
  //! \param[in] distance The distance between the particles.
  //! \param[in] id1 The id of the first particle.
  //! \param[in] id2 The id of the second particle.
  inline void LennardJones::forceStrength(RealType *normal, const RealType distance, const int id1, const int id2) const {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = Base::simData->type[id1]<0 ? 0 : 24.; //--
    c1 = Base::simData->type[id2]<0 ? 0 : c1; //--

    RealType *sg = simData->sg;
    // Calculate the magnitude
    RealType gamma = (sg[id1]+sg[id2]) / (distance*cutoff);
    RealType g3 = gamma*gamma*gamma, g6 = sqr(g3), g12 = sqr(g6);
    // The c1*c2 makes sure the magnitude is zero if either particles are of type -1
    RealType magnitude = c1*strength*(2.*g12 - g6)*(1./distance); 
    Force::virial += magnitude;
    // Make vectorial
    scalarMultVec(magnitude, normal);
    // Add forces
    plusEqVec (Base::simData->f[id1], normal);
    minusEqVec(Base::simData->f[id2], normal);
  }

}
