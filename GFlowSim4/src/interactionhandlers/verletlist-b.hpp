#ifndef __VERLET_LIST_B_HPP__GFLOW__
#define __VERLET_LIST_B_HPP__GFLOW__

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class VerletListB : public InteractionHandler {
    public:
    //! Constructor
    VerletListB(GFlow *gflow);

    //! @brief Add a pair of interacting particles.
    virtual void addPair(const int, const int);

    //! @brief Clear out all the data in the handler, allowing it to get new information.
    virtual void clear();

    //! @brief Returns a characteriztic size indicating the number of interactions the handler will test for. 
    //! 
    //! This may be the size of the array the particles are stored in, or the number of interactions the
    //! handler will try. If unknown, the function returns -1. This could occur if e.g. we are travering 
    //! through linked cells, and don't have an accurate guess of how many interactions there will be.
    virtual int size() const;

    //! @brief Return whether this handler needs to be constructed from the outside. 
    //!
    //! In other words, does "add pair" do anything. This is true by default.
    virtual bool needsConstruction() { return true; }

    //! @brief Iterate through interacting particles, executing the given kernel between them.
    //!
    //! @param kernel A function that is executed on all pairs of particles within cutoff distance
    //! of each other.
    //! @param param_pack Parameters used to evaluate the force.
    //! @param data_pack Data to be updated by the function.
    virtual void executeKernel(Kernel, const RealType*, RealType*) const;

  private:
    // Implementation

    //! @brief The particle ids
    vector<int> verlet;
    //! @brief The positions of the heads in the verlet list - this is mutable so we can terminate heads from the 
    //! (const) function executeKernel.
    mutable vector<int> heads;

    //! @brief The current head, for list creation
    int currentHead;

    //! @brief If true, we have added verlet.size() to heads, allowing us to use the heads list to 
    //! end our search through the verlet lists
    mutable bool terminated;
  };

}
#endif // __VERLET_LIST_B_HPP__GFLOW__