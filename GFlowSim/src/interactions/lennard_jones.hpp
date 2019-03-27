#ifndef __LENNARD_JONES_HPP__GFLOW__
#define __LENNARD_JONES_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class LennardJones : public Interaction {
  public:
    //! \brief Default constructor.
    LennardJones(GFlow *gflow, InteractionHandler *hndlr) : Interaction(gflow, hndlr), strength(DEFAULT_LENNARD_JONES_STRENGTH) { 
      cutoff = 2.5; 
    };

    //! \brief Set the cutoff factor.
    void setCutoff(RealType c) { cutoff = c>0 ? c : cutoff; }

    //! \brief Set the interaction strength.
    void setStrength(RealType s) { strength = s; }

  protected:

    //! \brief Strength of the LJ force.
    RealType strength;
  };

}
#endif // __LENNARD_JONES_HPP__GFLOW__