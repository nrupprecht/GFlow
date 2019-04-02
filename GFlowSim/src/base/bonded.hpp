#ifndef __BONDED_HPP__GFLOW__
#define __BONDED_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  /*
  *  \brief The base class for all bonded interactions.
  *
  */
  class Bonded : public Base {
  public:
    //! \brief Default constructor
    Bonded(GFlow*);

    //! \brief Calculate all the forces between bonded atoms.
    virtual void interact() const;

    // --- Accessors

    //! \brief Return the total number of bonds.
    virtual int size() const = 0;

    //! \brief Get the virial, used for calculating pressure
    RealType getVirial() const;

    //! \brief Get the potential energy.
    RealType getPotential() const;

    // --- Mutators

    //! \brief Set the do virial flag.
    void setDoVirial(bool);

    //! \brief Set the do potential flag.
    void setDoPotential(bool);

  private:
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
#endif // __BONDED_HPP__GFLOW__