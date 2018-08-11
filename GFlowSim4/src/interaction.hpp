#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include "gflow.hpp"
#include "interactionhandler.hpp"
#include "vectormath.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  /*
  *  \brief The base class for all interparticle forces.
  *
  *  A pair force between particles. This is the base class for forces. Forces keep a 
  *  verlet list of all the particles that might experience it.
  *
  */
  class Interaction : public Base {
  public:
    //! @brief Constructor
    Interaction(GFlow *);

    //! @brief Destructor
    ~Interaction();

    //! @brief Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() const;

    //! @brief Run an externally given kernel through this interaction's interaction handler
    virtual void executeKernel(Kernel, RealType*, RealType*) const;

    //! @brief Initialize for force calculation
    //!
    //! Resets the virial and returns whether forces should be calculated.
    bool initCalculation() const;

    // --- Accessors

    //! @brief Return the total length of the verlet list 
    int size() const;

    //! @brief Get the verlet list (get it as a const reference)
    const InteractionHandler* getInteractionHandler() const;

    //! @brief Get the virial, used for calculating pressure
    int getVirial() const;

    //! @brief Return whether this interaction needs to be constructed from the outside. 
    //!
    //! In other words, does "add pair" do anything. This is true by default.
    virtual bool needsConstruction() const;

    // --- Mutators

    //! @brief Clear this force's interaction handler
    void clear();

    //! @brief Add a pair pf particles - the first is the head
    void addPair(int, int);

    // GFlow is a friend class
    friend class GFlow;

  protected:

    //! @brief The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    InteractionHandler *handler;

    //! @brief A pointer to the force function
    Kernel forcePtr;

    //! @brief An array of force parameters. The length of this array will vary by force.
    RealType *parameters;

    //! @brief The virial, for calculating pressure.
    //!
    //! The pressure formula is: P = N k T/V + 1/(DIMENSIONS*V) \sum_i (r_i \dot F_i)
    //! This term should be used like: virial = \sum_i (r_i \dot F_i)
    mutable RealType virial;
  };

}
#endif // __FORCE_HPP__