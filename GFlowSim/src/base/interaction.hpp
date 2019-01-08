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
    //! \brief Constructor
    Interaction(GFlow *);

    //! \brief Destructor
    virtual ~Interaction();

    //! \brief Calculate all the forces between atoms in the verlet lists
    virtual void interact() const;

    //! \brief Compute and update interaction for pairs of particles.
    //!
    //! @param id1 Local id of the first particle.
    //! @param id2 Local id of the second particle.
    //! @param displacement The displacement from particle 1 to particle 2.
    //! @param distance The distance between the particles.
    virtual void compute(const int id1, const int id2, RealType *displacement, const RealType distance) const = 0;

    // --- Accessors

    //! \brief Return the total length of the verlet list.
    int size() const;

    //! \brief Get the verlet list (get it as a const reference)
    InteractionHandler* getInteractionHandler() const;

    //! \brief Get the virial, used for calculating pressure
    int getVirial() const;

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

    //! \brief An interaction kernel that can be used for the computation of particle interactions.
    Kernel kernel;

    //! \brief The virial, for calculating pressure.
    //!
    //! The pressure formula is: P = N k T/V + 1/(DIMENSIONS*V) \sum_i (r_i \dot F_i)
    //! This term should be used like: virial = \sum_i (r_i \dot F_i)
    mutable RealType virial = 0;

    //! \brief A pack of parameters.
    RealType *param_pack;
  };

}
#endif // __FORCE_HPP__
