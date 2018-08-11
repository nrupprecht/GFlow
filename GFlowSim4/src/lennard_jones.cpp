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
    // Add the force to the buffers
    plusEqVec (simData->f[id1], normal);
    minusEqVec(simData->f[id2], normal);
  }
  
}
