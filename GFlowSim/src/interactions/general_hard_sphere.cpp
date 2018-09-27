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
    int repulsionPlace = param_pack[0];
    int dissipationPlace = param_pack[1];
    int frictionPlace = param_pack[2];
    // Calculate force strength
    RealType repulsion = 0.5*(simData->dataF[repulsionPlace][id1] + simData->dataF[repulsionPlace][id2]);
    RealType magnitude = repulsion*(simData->sg[id1] + simData->sg[id2] - distance);
    // Update the virial with the normal component of force
    data_pack[0] += magnitude;
    // Normal force - force strength x Normal vector -> Sets normal to be the vectorial force between the particles.
    scalarMultVec(magnitude, normal);
    // Add the force to the buffers
    plusEqVec (simData->f[id1], normal);
    minusEqVec(simData->f[id2], normal);
  }

  void GeneralHardSphere::getPlaces() {
    parameters[0] = repulsionPlace   = Base::simData->getDataFEntry(repulsionLabel);
    parameters[1] = dissipationPlace = Base::simData->getDataFEntry(dissipationLabel);
    parameters[2] = frictionPlace    = Base::simData->getDataFEntry(frictionLabel);
  }

}