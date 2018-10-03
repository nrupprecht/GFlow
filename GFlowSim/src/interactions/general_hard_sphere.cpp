#include "general_hard_sphere.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  GeneralHardSphere::GeneralHardSphere(GFlow *gflow) : Interaction(gflow) {
    parameters = new RealType[3];
    // Set some default names 
    repulsionLabel   = "rp";
    dissipationLabel = "ds";
    frictionLabel    = "cf";
    // Get places
    getPlaces();
    // Set the force function
    kernelPtr = &force;
  };

  void GeneralHardSphere::initialize() {
    // Find data entries
    getPlaces();
    // Make sure there is angular data
    if (!Base::simData->usingAngularDynamics())
      throw DataEntryNotFound("GeneralHardSphere needs angular dynamics.");
    // Make sure all data entries were found - if not, throw an error.
    if (repulsionPlace<0 || dissipationPlace<0 || frictionPlace<0)
      throw DataEntryNotFound("GeneralHardSphere could not find a data entry.");
  }

  //! @param[in] normal
  //! @param[in] distance
  //! @param[in] id1
  //! @param[in] id2
  //! @param[in] simData
  //! @param[in] param_pack A parameter pack, passed in from force. Should be of the form { repulsion } (length 1).
  //! @param[in,out] data_pack Data to update in the function. Should be of the form  { virial } (length 1). 
  //! Add the f_i \dot r_i to this.
  void GeneralHardSphere::force(RealType* normal, const RealType distance, const int id1, const int id2, SimData *simData, 
    const RealType *param_pack, RealType *data_pack) 
  {
    int repulsionPlace   = param_pack[0];
    int dissipationPlace = param_pack[1];
    int frictionPlace    = param_pack[2];
    // Calculate force strength
    RealType repulsion = DEFAULT_HARD_SPHERE_REPULSION; // 0.5*(simData->dataF[repulsionPlace][id1] + simData->dataF[repulsionPlace][id2]);
    RealType Fn = repulsion*(simData->sg[id1] + simData->sg[id2] - distance);
    // Update the virial with the normal component of force
    data_pack[0] += Fn;
    
    // Frictional force - mu*Fn, damping force
    RealType mu = 0.1; // simData->dataF[frictionPlace][id1] * simData->dataF[frictionPlace][id2];
    RealType dissipation = 0.; // simData->dataF[dissipationPlace][id1] + simData->dataF[dissipationPlace][id2];
    // The direction of the relative tangential velocity gives the direction of the friction force
    RealType dV[DIMENSIONS];
    subtractVec(simData->v[id2], simData->v[id1], dV);
    // Normal velocity
    RealType Vn = dotVec(dV, normal);
    // Shear vector
    RealType shear[] = {normal[1], -normal[0]};
    // Shear velocity
    RealType Vs = dotVec(dV, shear) + simData->sg[id1]*simData->om[id1] + simData->sg[id2]*simData->om[id2];
    // Adjust magnitude of normal force for dissipation
    Fn += dissipation*clamp(Vn);
    RealType Fs = mu*Fn*sign(Vs);

    // Normal force - Fn x Normal vector -> Sets normal to be the normal force
    scalarMultVec(Fn, normal);
    // Shear force - Fs x shear vector -> Sets shear to be the shear force
    scalarMultVec(Fs, shear);
    // Total force - accumulate in normal
    plusEqVec(normal, shear);

    // Add the force to the buffers
    plusEqVec (simData->f[id1], normal);
    minusEqVec(simData->f[id2], normal);

    // Add torques
    simData->tq[id1] -= simData->sg[id1]*Fs;
    simData->tq[id2] -= simData->sg[id2]*Fs;
  }

  void GeneralHardSphere::getPlaces() {
    parameters[0] = repulsionPlace   = Base::simData->getDataFEntry(repulsionLabel);
    parameters[1] = dissipationPlace = Base::simData->getDataFEntry(dissipationLabel);
    parameters[2] = frictionPlace    = Base::simData->getDataFEntry(frictionLabel);
  }

}