#include "hard_sphere.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Interaction(gflow) {
    parameters = new RealType;
    *parameters = DEFAULT_HARD_SPHERE_REPULSION;
    // Set the force function
    kernelPtr = &force;
  };

  void HardSphere::setRepulsion(RealType r) { 
    parameters[0] = r; 
  }

  //! @param[in] normal
  //! @param[in] distance
  //! @param[in] id1
  //! @param[in] id2
  //! @param[in] simData
  //! @param[in] param_pack A parameter pack, passed in from force. Should be of the form { repulsion } (length 1).
  //! @param[in,out] data_pack Data to update in the function. Should be of the form  { virial } (length 1). 
  //! Add the f_i \dot r_i to this.
  void HardSphere::force(RealType* normal, const RealType distance, const int id1, const int id2, SimData *simData, 
    const RealType *param_pack, RealType *data_pack) 
  {
    // Calculate force strength
    RealType magnitude = param_pack[0]*(simData->sg[id1] + simData->sg[id2] - distance);
    // Update the virial
    data_pack[0] += magnitude;
    // Force strength x Normal vector -> Sets normal to be the vectorial force between the particles.
    scalarMultVec(magnitude, normal);
    // Add the force to the buffers
    plusEqVec (simData->f[id1], normal);
    minusEqVec(simData->f[id2], normal);
  }

}
