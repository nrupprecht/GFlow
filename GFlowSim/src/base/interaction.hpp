#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include "../gflow.hpp"
#include "interactionhandler.hpp"
#include "../utility/vectormath.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  template<typename float_type> struct Particle_Pack {
    float_type *X1, *X2;
    float_type *V1, *V2;
    float_type *F1, *F2;
    float_type *normals;
    //! @brief Cutoff mask is 0 or 1 based on whether the particles are within cutoff distance
    float_type cutoff_mask; 

    //! @brief Extra data
    float_type *dataF;
  };

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

    //! @brief Initialize the force, check if all the special data (dataF, dataI) the force needs exists, make
    //! sure parameter packs are up to date.
    virtual void initialize()=0;

    //! @brief Calculate all the forces between atoms in the verlet lists
    virtual void interact() const;

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
    //! In other words, does "add pair" do anything. It is actually determined by the interaction handler.
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
    Kernel kernelPtr;

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