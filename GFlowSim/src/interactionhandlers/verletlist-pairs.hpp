#ifndef __VERLET_LIST_PAIRS__
#define __VERLET_LIST_PAIRS__

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class VerletListPairs : public InteractionHandler {
  public:
    VerletListPairs(GFlow*);

    // --- Inherited members

    //! \brief Add a pair of interacting particles.
    virtual void addPair(const int, const int) override;

    //! \brief This function lets the handler know that no more pairs will be added to it. 
    //! For VerletListPairs, this function does nothing.
    virtual void close() override {};

    //! \brief Clear the vectors of data.
    virtual void clear() override;

    //! \brief Return the total length of the verlet list.
    virtual int size() const override;

    //! \brief The verlet list
    vector<int> verlet;

  };

}
#endif // __VERLET_LIST_PAIRS__