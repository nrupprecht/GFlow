#include "coulombforce.hpp"

namespace GFlowSimulation {

  CoulumbForce::CoulumbForce(GFlow *gflow) : Interaction(gflow) {
    // Set up parameters
    parameters = new RealType[2];
    parameters[0] = 1.;
    parameters[1] = 2.5;
    // Set the force function
    kernelPtr = force;
  }

  void CoulumbForce::force(RealType* normal, const RealType distance, const int id1, const int id2, SimData *simData, 
    const RealType *param_pack, RealType *data_pack)
  {
    // Calculate the magnitude of the force
    RealType magnitude = param_pack[0]*sqr((simData->sg[id1] + simData->sg[id2])/(param_pack[1]*distance));
    // Update virial
    data_pack[0] += magnitude;
    // Make vectorial
    scalarMultVec(magnitude, normal);
    // Add the force to the buffers
    plusEqVec (simData->f[id1], normal);
    minusEqVec(simData->f[id2], normal);
  }

}