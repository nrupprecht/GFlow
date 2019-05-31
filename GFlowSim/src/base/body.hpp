#ifndef __BODY_HPP__GFLOW__
#define __BODY_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class Body : public Base {
  public:
    //! \brief Default constructor.
    Body(GFlow *gflow) : Base(gflow) {};

    //! \brief Do whatever corrections the body particles need.
    virtual void correct() = 0;

    //! \brief Get the virial, used for calculating pressure
    RealType getVirial() const {
      return virial;
    }

    //! \brief Get the potential energy.
    RealType getPotential() const {
      return potential;
    }

    // --- Mutators

    //! \brief Set the do virial flag.
    void setDoVirial(bool d) {
      do_virial = d;
    }

    //! \brief Set the do potential flag.
    void setDoPotential(bool d) {
      do_potential = d;
    }

  protected:

    //! \brief The virial, for calculating pressure.
    //!
    //! The pressure formula is: P = N k T/V + 1/(DIMENSIONS*V) \sum_i (r_i \dot F_i)
    //! This term should be used like: virial = \sum_i (r_i \dot F_i)
    mutable RealType virial = 0;

    //! \brief The potential energy for the interaction.
    mutable RealType potential = 0;

    //! \brief Whether to do virial calculation.
    bool do_virial = true;

    //! \brief Whether to do potential energy calculation.
    bool do_potential = true;
  };

}
#endif // __BODY_HPP__GFLOW__