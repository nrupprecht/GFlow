#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include "../gflow.hpp"
#include "interactionhandler.hpp"
#include "../utility/vectormath.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  /*
  *  \brief The base class for all interparticle forces and other interactions.
  *
  *  A pair force between particles. This is the base class for forces. Forces have an interaction handler that
  *  is responsible for storing the particles that the interaction should act on.
  *
  */
  class Interaction : public Base {
  public:
    //! \brief Default Constructor.
    Interaction(GFlow*);

    //! \brief Calculate all the forces between atoms in the verlet lists.
    virtual void interact() const;

    // --- Accessors

    //! \brief Get the cutoff for the interaction.
    RealType getCutoff() const;

    //! \brief Get the virial, used for calculating pressure
    RealType getVirial() const;

    //! \brief Get the potential energy.
    RealType getPotential() const;

    // --- Mutators

    //! \brief Set the do virial flag.
    void setDoVirial(bool);

    //! \brief Set the do potential flag.
    void setDoPotential(bool);

    //! \brief Add a pair of interacting particles.
    virtual void addPair(const int, const int);

    //! \brief Signals that the pair additions are done.
    virtual void close() {};

    //! \brief Clears the verlet list.
    virtual void clear();

    //! \brief Returns the size of the verlet list.
    virtual int size() const;

    // GFlow is a friend class
    friend class GFlow;

  protected:

    //! \brief The verlet list
    vector<int> verlet;

    //! \brief The cutoff, i.e. how many times greater than the particle radius can a particle interact with other particles.
    //!
    //! For example, for hard sphere interactions, the cutoff is 1. For Lennard jones particles, the cutoff is generally 2.5 or so.
    //! The default value is 1.
    RealType cutoff = 1.;

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
#endif // __FORCE_HPP__
