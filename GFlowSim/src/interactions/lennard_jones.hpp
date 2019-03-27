#ifndef __LENNARD_JONES_HPP__GFLOW__
#define __LENNARD_JONES_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class LennardJones : public Interaction {
  public:
    //! \brief Default constructor.
    LennardJones(GFlow *gflow, InteractionHandler *hndlr) : Interaction(gflow, hndlr), strength(DEFAULT_LENNARD_JONES_STRENGTH) { 
      setCutoff(2.5);
    };

    //! \brief Set the cutoff factor.
    void setCutoff(RealType c) { 
      cutoff = c>0 ? c : cutoff; 
      // Calculate the potential when r is the cutoff radius
      RealType gamma = 1./cutoff;
      RealType g3  = gamma*gamma*gamma; 
      RealType g6  = g3*g3;
      RealType g12 = g6*g6;
      // Set the potential energy shift
      potential_energy_shift = 4.*strength*(g12 - g6);
    }

    //! \brief Set the interaction strength.
    void setStrength(RealType s) { strength = s; }

  protected:

    //! \brief Strength of the LJ force.
    RealType strength;

    //! \brief Energy shift caused by cutoff.
    RealType potential_energy_shift = 0;
  };

}
#endif // __LENNARD_JONES_HPP__GFLOW__