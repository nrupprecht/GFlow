#ifndef __LENNARD_JONES__GFLOW__
#define __LENNARD_JONES__GFLOW__

#include "interaction.hpp"

namespace GFlowSimulation {

  /*
  *  \brief LennardJones where all particles have the same force strength.
  *
  *  Lennard Jones force. The particle sigma will represent the force cutoff,
  *  generally 2.5*sig, where sig is the inter-particle distance where V=0.
  *  We use cutoff=2.5 by default, but it can be changed. Strength is the
  *  "epsilon" parameter in LJ.
  *
  */
  class LennardJones : public Interaction {
  public:
    //! @brief Constructor
    LennardJones(GFlow *);

    //! @brief Set the lennard jones interaction strength. Must be non-negative.
    void setStrength(RealType);

    //! @brief Set the lennard jones cutoff range. Must be at least 1.
    void setCutoff(RealType);

    //! @param[in] normal
    //! @param[in] distance
    //! @param[in] id1
    //! @param[in] id2
    //! @param[in] simData
    //! @param[in] param_pack A parameter pack, passed in from force. Contains characteristic 
    //! constants of the force, and extra data the force needs.
    //! @param[in,out] data_pack Data to be updated by the function.
    static void force(RealType*, const RealType, const int, const int, const SimData*, const RealType*, RealType*);

  private:
    //! @brief Calculate force strength
    void forceStrength(RealType*, const RealType, const int, const int) const;

    // The parameters for LJ are the LJ strength (parameters[0]), and the cuttoff (parameters[1]).
  };

}
#endif // __LENNARD_JONES__GFLOW__