#include "lennard_jones.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Force(gflow) {
    parameters = new RealType[2];
    parameters[0] = DEFAULT_LENNARD_JONES_STRENGTH;
    parameters[1] = DEFAULT_LENNARD_JONES_CUTOFF;
    // Set the force function
    forcePtr = force;
  };

  /*
  void LennardJones::calculateForces() const {
    // Check if there are forces to calculate. Virial is reset.
    if (!Force::initCalculation()) return; 

    //RealType param_pack[] = { strength, cutoff };
    //verletList->forceKernel(&force, parameters, &virial);

    // Get the data we need
    int nverlet = verletList->vlSize(), id1, id2; // List length, id pointers
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // Get verlet list data
    const int *verlet = verletList->getVerlet();

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
  */

  void LennardJones::setStrength(RealType s) {
    if (s>=0) parameters[0] = s;
  }

  void LennardJones::setCutoff(RealType cut) {
    if (cut>1.) parameters[1] = cut;
  }

  //! @param[in,out] normal
  //! @param[in] distance
  //! @param[in] id1
  //! @param[in] id2
  //! @param[in] simData
  //! @param[in] param_pack A parameter pack, passed in from force. Should be of the form { strength, cutoff } (length 2).
   //! @param[in,out] data_pack Data to update in the function. Should be of the form  { virial } (length 1). 
  //! Add the f_i \dot r_i to this.
  void LennardJones::force(RealType* normal, const RealType distance, const int id1, const int id2, const SimData *simData, 
    const RealType *param_pack, RealType *data_pack)
  {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = (simData->type[id1]<0 || simData->type[id2]<0) ? 0 : 24.; //--

    // Calculate the magnitude
    RealType gamma = (simData->sg[id1]+simData->sg[id2]) / (distance*param_pack[1]);
    RealType g3 = gamma*gamma*gamma, g6 = sqr(g3), g12 = sqr(g6);
    // Debugging assertion
    #if DEBUG==1
    assert(distance>0)
    #endif
    // The c1 factor makes sure the magnitude is zero if either particles are of type -1
    RealType magnitude = c1*param_pack[0]*(2.*g12 - g6)*(1./distance); 
    data_pack[0] += magnitude;
    // Make vectorial
    scalarMultVec(magnitude, normal);
  }

  //! \param[in] normal The normal of the displacement between the particles. Points from id2 to id1 (?). It is safe
  //! to change this parameter.
  //! \param[in] distance The distance between the particles.
  //! \param[in] id1 The id of the first particle.
  //! \param[in] id2 The id of the second particle.
  inline void LennardJones::forceStrength(RealType *normal, const RealType distance, const int id1, const int id2) const {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = (simData->type[id1]<0 || simData->type[id2]<0) ? 0 : 24.; //--

    // Calculate the magnitude
    RealType gamma = (Base::simData->sg[id1]+Base::simData->sg[id2]) / (distance*parameters[1]);
    RealType g3 = gamma*gamma*gamma, g6 = sqr(g3), g12 = sqr(g6);
    // Debugging assertion
    #if DEBUG==1
    assert(distance>0)
    #endif
    // The c1 factor makes sure the magnitude is zero if either particles are of type -1
    RealType magnitude = c1*parameters[0]*(2.*g12 - g6)*(1./distance); 
    Force::virial += magnitude;
    // Make vectorial
    scalarMultVec(magnitude, normal);
    // Add forces
    plusEqVec (Base::simData->f[id1], normal);
    minusEqVec(Base::simData->f[id2], normal);
  }

}
