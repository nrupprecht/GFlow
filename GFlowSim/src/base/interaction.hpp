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
    //! @brief Constructor
    Interaction(GFlow *);

    //! @brief Destructor
    ~Interaction();

    //! @brief Initialize the force, check if all the special data (dataF, dataI) the force needs exists, make
    //! sure parameter packs are up to date.
    virtual void initialize()=0;

    //! @brief Calculate all the forces between atoms in the verlet lists
    virtual void interact() const;

    //! @brief Run an externally given kernel through this interaction's interaction handler
    virtual void executeKernel(Kernel<simd_float>, Kernel<float>, const RealType*, RealType*, const vector<int>&) const;

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
    //! In other words, does "add pair" do anything. It is actually determined by the interaction handler.
    virtual bool needsConstruction() const;

    // --- Mutators

    //! @brief Clear this force's interaction handler
    void clear();

    //! @brief Add a pair pf particles - the first is the head
    virtual void addPair(int, int);

    //! @brief Signals that the pair additions are done.
    virtual void close();

    // GFlow is a friend class
    friend class GFlow;

  protected:

    //! @brief The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    InteractionHandler *handler;

    //! @brief A pointer to the force function
    Kernel<simd_float> simd_kernelPtr = nullptr;
    Kernel<float> serial_kernelPtr = nullptr;

    //! @brief An array of force parameters. The length of this array will vary by force.
    RealType *parameters = nullptr;

    //! @brief The pointers to the arrays of data in simdata that should be packed for the interaction kernel function.
    vector<int> data_needed;

    //! @brief The virial, for calculating pressure.
    //!
    //! The pressure formula is: P = N k T/V + 1/(DIMENSIONS*V) \sum_i (r_i \dot F_i)
    //! This term should be used like: virial = \sum_i (r_i \dot F_i)
    mutable RealType virial = 0;
  };

}
#endif // __FORCE_HPP__