#ifndef __LENNARD_JONES_HPP__GFLOW__
#define __LENNARD_JONES_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class LennardJones : public Interaction {
  public:
    //! \brief Default constructor.
    LennardJones(GFlow *gflow, InteractionHandler *hndlr) 
      : Interaction(gflow, hndlr), cutoff(DEFAULT_LENNARD_JONES_CUTOFF), strength(DEFAULT_LENNARD_JONES_STRENGTH) {};

    //! \brief Set the cutoff factor.
    void setCutoff(RealType c) { cutoff = c>0 ? c : cutoff; }

    //! \brief Set the interaction strength.
    void setStrength(RealType s) { strength = s; }

  protected:  
    //! \brief The cutoff factor. Generally, this is 2.5 times the "radius" of the particle.
    RealType cutoff;

    //! \brief Strength of the LJ force.
    RealType strength;
  };

}
#endif // __LENNARD_JONES_HPP__GFLOW__