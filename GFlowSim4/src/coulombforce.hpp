#ifndef __COULOMB_FORCE_HPP__GFLOW__
#define __COULOMB_FORCE_HPP__GFLOW__

#include "interaction.hpp"

namespace GFlowSimulation {

  /** 
  *  @brief A coulombic force where all particles have the same charge and force contant.
  *
  *  The force between particles is proportional to 1/r^2, with a proportionality 
  *  constant stored by the force. Charges are all taken to be 1.
  *
  *  The (charge^2 * strength) is parameters[0]
  */
  class CoulumbForce : public Interaction {
  public:
    //! Constructor
    CoulumbForce(GFlow*);

    //! Set the repulsion parameter.
    void setConstant(RealType);

    //! @param[in] normal
    //! @param[in] distance
    //! @param[in] id1
    //! @param[in] id2
    //! @param[in] simData
    //! @param[in] param_pack A parameter pack, passed in from force. Contains characteristic 
    //! constants of the force, and extra data the force needs. Of the form { strength, cutoff } (length 2).
    //! @param[in,out] data_pack Data to be updated by the function.
    static void force(RealType*, const RealType, const int, const int, SimData*, const RealType*, RealType*);
  };
  
}
#endif // __COULOMB_FORCE_HPP__GFLOW__