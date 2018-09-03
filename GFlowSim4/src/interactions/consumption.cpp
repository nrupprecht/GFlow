#include "consumption.hpp"

namespace GFlowSimulation {

  Consumption::Consumption(GFlow *gflow) : Interaction(gflow) {
    //parameters = new RealType;
    //*parameters = DEFAULT_HARD_SPHERE_REPULSION;
    
    // Set the force function
    kernelPtr = &consume;
  }

  //! @param[in] normal
  //! @param[in] distance
  //! @param[in] id1
  //! @param[in] id2
  //! @param[in] simData
  //! @param[in] param_pack A parameter pack, passed in from force. Should be of the form { repulsion } (length 1).
  //! @param[in,out] data_pack Data to update in the function. Should be of the form  { virial } (length 1). 
  //! Add the f_i \dot r_i to this.
  void Consumption::consume(RealType* normal, const RealType distance, const int id1, const int id2, SimData *simData, 
    const RealType *param_pack, RealType *data_pack) 
  {
    RealType random = drand48();
    const RealType probability = 0.001;
    // Smaller type eats larger type
    if (simData->type[id1]<simData->type[id2] && random<probability)
      simData->markForRemoval(id2);
    else if (random<probability) 
      simData->markForRemoval(id1);
  }

}