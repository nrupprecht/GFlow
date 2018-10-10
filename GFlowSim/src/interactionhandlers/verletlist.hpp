#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class VerletList : public InteractionHandler {
  public:
    //! @brief Constructor.
    VerletList(GFlow*);

    // --- Inherited members

    //! @brief Add a pair of interacting particles.
    virtual void addPair(const int, const int) override;

    //! @brief Add a last entry to the heads list to denote the end of the verlet list (a sentinal value).
    virtual void close() override;

    //! @brief Clear the vectors of data.
    virtual void clear();

    //! @brief Return the total length of the verlet list.
    virtual int size() const;

    //! @brief Iterate through interacting particles, executing the given kernel between them.
    //!
    //! @param kernel A function that is executed on all pairs of particles within cutoff distance
    //! of each other.
    //! @param param_pack Parameters used to evaluate the force.
    //! @param data_pack Data to be updated by the function.
    virtual void executeKernel(Kernel, const RealType*, RealType*, const vector<int>&) const override;

  private:

    //! @brief The verlet list
    vector<int> verlet;
    vector<int> heads;
    int lastHead;

    //! @brief The minimum number of particles in a list for which we will use simd instead of serial
    int min_simd_size;

    //! @brief If true, we can use simd functions
    bool use_simd;
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__
