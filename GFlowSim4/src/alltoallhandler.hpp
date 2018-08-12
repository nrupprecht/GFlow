#ifndef __ALL_TO_ALL_INTERACTION_HPP__GFLOW__
#define __ALL_TO_ALL_INTERACTION_HPP__GFLOW__

#include "interactionhandler.hpp"

namespace GFlowSimulation {

  class AllToAllHandler : public InteractionHandler {
  public:
    //! @brief Constructor.
    AllToAllHandler(GFlow*);

    //! @brief This function does nothing, though we must implement it.
    virtual void addPair(const int, const int) {};

    //! @brief This function does nothing, though we must implement it.
    virtual void clear() {};

    //! @brief We will check all particle pairs for interactions
    virtual int size() const;

    //! @brief Return whether this handler needs to be constructed from the outside. 
    //!
    //! In other words, does "add pair" do anything. This is false, since we are just going to 
    //! loop over all particle pairs
    virtual bool needsConstruction() { return false; }

    //! @brief Iterate through interacting particles, executing the given kernel between them.
    //!
    //! @param kernel A function that is executed on all pairs of particles within cutoff distance
    //! of each other.
    //! @param param_pack Parameters used to evaluate the force.
    //! @param data_pack Data to be updated by the function.
    virtual void executeKernel(Kernel, const RealType*, RealType*) const;

  };

}
#endif // __ALL_TO_ALL_INTERACTION_HPP__GFLOW__