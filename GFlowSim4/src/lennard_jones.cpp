#include "lennard_jones.hpp"

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Interaction(gflow) {
    // Set up parameters
    parameters = new RealType[2];
    parameters[0] = DEFAULT_LENNARD_JONES_STRENGTH;
    parameters[1] = DEFAULT_LENNARD_JONES_CUTOFF;
    // Set the force function
    kernelPtr = force;
  };

  void LennardJones::setStrength(RealType s) {
    if (s>=0) parameters[0] = s;
  }

  void LennardJones::setCutoff(RealType cut) {
    if (cut>1.) parameters[1] = cut;
  }

  // param_pack Should be of the form { strength, cutoff } (length 2).
  // data_pack Should be of the form  { virial } (length 1). 
  void LennardJones::force(RealType* normal, const RealType distance, const int id1, const int id2, SimData *simData, 
    const RealType *param_pack, RealType *data_pack)
  {
    // Calculate the magnitude of the force
    RealType gamma = (simData->sg[id1]+simData->sg[id2]) / (distance*param_pack[1]);
    RealType g3 = gamma*gamma*gamma, g6 = sqr(g3), g12 = sqr(g6);
    RealType magnitude = 24*param_pack[0]*(2.*g12 - g6)*(1./distance); 
    // Update the virial
    data_pack[0] += magnitude;
    // Make vectorial
    scalarMultVec(magnitude, normal);
    // Add the force to the buffers
    plusEqVec (simData->f[id1], normal);
    minusEqVec(simData->f[id2], normal);
  }
  
}
