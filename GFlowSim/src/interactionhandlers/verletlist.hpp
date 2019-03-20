#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  class VerletList : public InteractionHandler {
  public:
    //! \brief Constructor.
    VerletList(GFlow*);

    // --- Inherited members

    //! \brief Add a pair of interacting particles.
    virtual void addPair(const int, const int) override;

    //! \brief Add a last entry to the heads list to denote the end of the verlet list (a sentinal value).
    virtual void close() override;

    //! \brief Clear the vectors of data.
    virtual void clear() override;

    //! \brief Return the total length of the verlet list.
    virtual int size() const override;

    //! \brief The verlet list.
    vector<int> verlet;
    
    //! \brief Pointers to heads in the verlet list.
    vector<int> heads;

  private:
    //! \brief The current head in the verlet list.
    int lastHead;
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__
