#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include "gflow.hpp"
#include "interactionhandler.hpp"

namespace GFlowSimulation {

  /*
  *  \brief The base class for all interparticle forces.
  *
  *  A pair force between particles. This is the base class for forces. Forces keep a 
  *  verlet list of all the particles that might experience it.
  *
  */
  class Force : public Base {
  public:
    //! Constructor
    Force(GFlow *);

    //! Destructor
    ~Force();

    //! Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() const;

    //! Run an externally given kernel through this interaction's interaction handler
    virtual void executeKernel(Kernel, RealType*, RealType*) const;

    //! \brief Initialize for force calculation
    //!
    //! Resets the virial and returns whether forces should be calculated.
    bool initCalculation() const;

    // --- Accessors

    //! Return the total length of the verlet list 
    int size() const;

    //! Get the verlet list (get it as a const reference)
    const InteractionHandler* getInteractionHandler() const;

    //! Get the virial, used for calculating pressure
    int getVirial() const;

    // --- Mutators

    //! Clear this force's interaction handler
    void clear();

    //! Add a pair pf particles - the first is the head
    void addPair(int, int);

    // GFlow is a friend class
    friend class GFlow;

  protected:

    //! The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    InteractionHandler *handler;

    //! A pointer to the force function
    Kernel forcePtr;

    //! An array of force parameters. The length of this array will vary by force.
    RealType *parameters;

    //! The virial, for calculating pressure.
    //! The pressure formula is: P = N k T/V + 1/(DIMENSIONS*V) \sum_i (r_i \dot F_i)
    //! This term should be used like: virial = \sum_i (r_i \dot F_i)
    mutable RealType virial;
  };

}
#endif // __FORCE_HPP__