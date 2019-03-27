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
  *  A pair force between particles. This is the base class for forces. Forces keep a 
  *  verlet list of all the particles that might experience it.
  *
  */
  class Interaction : public Base {
  public:
    //! \brief Default Constructor.
    Interaction(GFlow *);

    //! \brief Constructor that sets the interaction handler.
    Interaction(GFlow *, InteractionHandler*);

    //! \brief Destructor.
    virtual ~Interaction();

    //! \brief Calculate all the forces between atoms in the verlet lists.
    virtual void interact() const;

    // --- Accessors

    //! \brief Return the total length of the verlet list.
    int size() const;

    //! \brief Get the verlet list (get it as a const reference)
    InteractionHandler* getInteractionHandler() const;

    //! \brief Get the cutoff for the interaction.
    RealType getCutoff() const;

    //! \brief Get the virial, used for calculating pressure
    RealType getVirial() const;

    //! \brief Get the potential energy.
    RealType getPotential() const;

    // --- Mutators

    //! \brief Clear this force's interaction handler
    void clear();

    //! \brief Add a pair pf particles - the first is the head
    virtual void addPair(int, int);

    //! \brief Signals that the pair additions are done.
    virtual void close();

    // GFlow is a friend class
    friend class GFlow;

  protected:

    //! \brief The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    InteractionHandler *handler;

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
