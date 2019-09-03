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
  *  A pair force between particles. This is the base class for forces.
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

    //! \brief Suggests a safe timescale given the minimum mass of a particle that has this interaction.
    //!
    //! Corresponds to e.g. the period of an interaction (for spring forces). 
    //! Returns -1 if no guidance is given.
    virtual RealType suggest_timescale(RealType) const;

    // --- Mutators

    //! \brief Set the do virial flag.
    void setDoVirial(bool);

    //! \brief Set the do potential flag.
    void setDoPotential(bool);

    //! \brief Add a pair of particles whose distance may need to be wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPair(const int, const int, const int);

    //! \brief Signals that the pair additions are done.
    virtual void close() {};

    //! \brief Clears the verlet list.
    virtual void clear();

    //! \brief Returns the number of interaction pairs in the interaction.
    virtual int size() const;

    // GFlow is a friend class
    friend class GFlow;

  protected:

    //! \brief The verlet lists for interactions that do not need to worry about wrapping distances.
    //!
    //! There are three types of verlet lists: \n
    //!   0 - Interactions that don't need the minimum image convention.
    //!   1 - Interactions that do need the minimum image convention.
    //!   2 - Interaction list for ghost particles.
    //!
    //! Perhaps one day, this will be extended so there can be more verlet lists, but for now, this should be fine.
    vector<int> verlet[3];

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
