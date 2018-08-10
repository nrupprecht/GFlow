#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

#include "gflow.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  typedef void (*ForceKernel) (RealType*, const RealType, const int, const int, const class SimData*, const RealType*, RealType*);

  class VerletList : public Base {
  public:
    //! @brief Constructor.
    VerletList(GFlow*);

    //! @brief Copy Constructor.
    VerletList(const VerletList&);

    //! @brief Destructor.
    ~VerletList();

    //! @brief Equals operator.
    VerletList& operator=(const VerletList&);

    // --- Mutators

    //! @brief Add a pair of interacting particles.
    void addPair(const int, const int);

    //! @brief Set sizes (but not capacities) to zero, effectively "clearing" out the data.
    void clear();

    //! @brief Return the total length of the verlet list.
    int vlSize() const;

    //! @brief Get a (const) pointer to the verlet array.
    const int* getVerlet() const;

    //! @brief THIS IS A TEST
    //! @param force A function static function that tells how to evaluate the force between particles.
    //! @param param_pack Parameters used to evaluate the force.
    //! @param data_pack Data to be updated by the function.
    void forceKernel(ForceKernel, const RealType*, RealType*) const;

  private:

    // --- Helper functions
    inline void resizeVerlet(); 

    // --- Data
    int *verlet;
    int vsize, vcapacity;
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__
