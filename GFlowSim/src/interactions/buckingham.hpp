#ifndef __BUCKINGHAM_HPP__GFLOW__
#define __BUCKINGHAM_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class Buckingham : public Interaction {
  public:
    //! \brief Default constructor.
    Buckingham(GFlow *gflow, InteractionHandler *hndlr) 
      : Interaction(gflow, hndlr), ratio(1.), strength(DEFAULT_LENNARD_JONES_STRENGTH) {
      cutoff = 2.5; 
    };

    //! \brief Ratio setting constructor.
    Buckingham(GFlow *gflow, RealType r, InteractionHandler *hndlr) 
      : Interaction(gflow, hndlr), ratio(r), strength(DEFAULT_LENNARD_JONES_STRENGTH) {
      cutoff = 2.5;   
    };

    //! \brief Set the r1 parameter.
    void setRatio(RealType r) { ratio = r; }

    //! \brief Set the strength parameter.
    void setStrength(RealType s) { strength = s; }

  protected:
    //! \brief The ratio r0/r1 factor, where r0 is the sum of the particle radii.
    RealType ratio;

    //! \brief The strength of the interaction.
    RealType strength;

    //! \brief The minimum allowable force (a large negative force). The force that occurs at the inner cutoff
    RealType inner_force = -1.;
  };

}
#endif // __BUCKINGHAM_HPP__GFLOW__