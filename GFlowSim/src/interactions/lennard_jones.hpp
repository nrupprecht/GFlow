#ifndef __LENNARD_JONES__GFLOW__
#define __LENNARD_JONES__GFLOW__

#include "../base/interaction.hpp"
#include "../utility/simd_generic.hpp"

namespace GFlowSimulation {

  /**
  *  \brief LennardJones where all particles have the same force strength.
  *
  *  Lennard Jones force. The particle sigma will represent the force cutoff,
  *  generally 2.5*sig, where sig is the inter-particle distance where V=0.
  *  We use cutoff=2.5 by default, but it can be changed. Strength is the
  *  "epsilon" parameter in LJ.
  *
  *  The parameters for LJ are the LJ strength (parameters[0]), and the cuttoff (parameters[1]).
  */
  class LennardJones : public Interaction {
  public:
    //! @brief Constructor
    LennardJones(GFlow *);

    //! @brief Set the lennard jones interaction strength. Must be non-negative.
    void setStrength(RealType);

    //! @brief Set the lennard jones cutoff range. Must be at least 1.
    void setCutoff(RealType);

    virtual void compute(const int, const int, RealType*, const RealType) const override;

    //! \brief An interaction kernel.
    static void kernel(SimData*, int, int, RealType*, RealType, RealType*, int);

  private:
    //! @brief The LJ force strength.
    RealType strength;
    //! @brief The LJ force cutoff.
    RealType cutoff;
  };

}
#endif // __LENNARD_JONES__GFLOW__