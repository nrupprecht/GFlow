#ifndef __INTERACTION_HANDLER_HPP__GFLOW__
#define __INTERACTION_HANDLER_HPP__GFLOW__

#include "../gflow.hpp"
#include "../utility/vectormath.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Sets up interactions.
  *
  *  Can, for example, create linked cells or verlet lists for helping interactions.
  */
  class InteractionHandler : public Base {
  public:
    //! Constructor
    InteractionHandler(GFlow *gflow) : Base(gflow) {};

    //! \brief Construct objects for interactions.
    virtual void construct();

  private:

    //! \brief Interaction grid.
    vector<vector<Interaction*> > grid;

    //! \brief Cutoff factors for each particle type
    vector<RealType> max_cutoffs;

  };

}
#endif // __INTERACTION_HANDLER_HPP__GFLOW__
