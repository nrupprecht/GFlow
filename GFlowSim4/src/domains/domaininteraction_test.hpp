#ifndef __DOMAIN_INTERACTION_TEST_HPP__GFLOW__
#define __DOMAIN_INTERACTION_TEST_HPP__GFLOW__

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class DomainInteraction : public InteractionHandler {
  public:
    //! @brief Constructor.
    DomainInteraction(GFlow*);

    //! @brief Add a pair of interacting particles.
    virtual void addPair(const int, const int) {};

    //! @brief Set sizes (but not capacities) to zero, effectively "clearing" out the data.
    virtual void clear() {};

    //! @brief Return the total length of the verlet list.
    virtual int size() const { return -1; };

    //! @brief Return whether this handler needs to be constructed from the outside. 
    //!
    //! In other words, does "add pair" do anything. This is true by default.
    virtual bool needsConstruction() { return false; }

    //! @brief Iterate through interacting particles, executing the given kernel between them.
    //!
    //! @param kernel A function that is executed on all pairs of particles within cutoff distance
    //! of each other.
    //! @param param_pack Parameters used to evaluate the force.
    //! @param data_pack Data to be updated by the function.
    virtual void executeKernel(Kernel, const RealType*, RealType*) const;

  private:

  };

}
#endif // __DOMAIN_INTERACTION_TEST_HPP__GFLOW__