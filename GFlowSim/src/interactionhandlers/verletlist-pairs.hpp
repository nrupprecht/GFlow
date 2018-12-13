#ifndef __VERLET_LIST_PAIRS__
#define __VERLET_LIST_PAIRS__

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class VerletListPairs : public InteractionHandler {
  public:
    VerletListPairs(GFlow*);

    // --- Inherited members

    //! @brief Add a pair of interacting particles.
    virtual void addPair(const int, const int) override;

    //! @brief Add a last entry to the heads list to denote the end of the verlet list (a sentinal value).
    virtual void close() override {};

    //! @brief Clear the vectors of data.
    virtual void clear() override;

    //! @brief Return the total length of the verlet list.
    virtual int size() const override;

    //! @brief Iterate through the pairs of interacting particles, hand them to Interaction's compute function.
    virtual void execute(const Interaction*) const override;

  private:
    //! @brief The verlet list
    vector<int> verlet;

  };

}
#endif // __VERLET_LIST_PAIRS__