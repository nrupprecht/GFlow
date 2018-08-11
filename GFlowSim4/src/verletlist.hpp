#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

#include "interactionhandler.hpp"

namespace GFlowSimulation {

  class VerletList : public InteractionHandler {
  public:
    //! @brief Constructor.
    VerletList(GFlow*);

    //! @brief Copy Constructor.
    VerletList(const VerletList&);

    //! @brief Destructor.
    ~VerletList();

    //! @brief Equals operator.
    VerletList& operator=(const VerletList&);

    // --- Inherited members

    //! @brief Add a pair of interacting particles.
    virtual void addPair(const int, const int);

    //! @brief Set sizes (but not capacities) to zero, effectively "clearing" out the data.
    virtual void clear();

    //! @brief Return the total length of the verlet list.
    virtual int size() const;

    //! @brief Iterate through interacting particles, executing the given kernel between them.
    //!
    //! @param kernel A function that is executed on all pairs of particles within cutoff distance
    //! of each other.
    //! @param param_pack Parameters used to evaluate the force.
    //! @param data_pack Data to be updated by the function.
    virtual void executeKernel(Kernel, const RealType*, RealType*) const;

    // --- Accessors

    //! @brief Get a (const) pointer to the verlet array.
    const int* getVerlet() const;

  private:

    // --- Helper functions
    inline void resizeVerlet(); 

    // --- Data
    int *verlet;
    int vsize, vcapacity;
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__
