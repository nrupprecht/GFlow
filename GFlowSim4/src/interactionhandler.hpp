#ifndef __INTERACTION_HANDLER_HPP__GFLOW__
#define __INTERACTION_HANDLER_HPP__GFLOW__

#include "gflow.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  //! @brief This defines what type of force kernel function the interaction handler expects to be passed in to it.
  typedef void (*ForceKernel) (RealType*, const RealType, const int, const int, const class SimData*, const RealType*, RealType*);


  /**
  *  @brief Handles the storage and traversal of interacting particles.
  *
  *  This is a generalization of what was the original verletlist class. In a more
  *  general setting, we may want a way of handling the particles that is not just
  *  verlet lists, e.g. linked cell based force handling. 
  *
  *  Children of this (purely abstract) class store the data on how particles interact
  *  however they choose, and know how to iterate through that data. The force passes
  *  in a force kernal that defines the force between pairs of particles. The 
  *  InteractionHandler goes through all the particles, applying the force kernal to 
  *  all particles i,j within sg[i] + sg[j] of one another.
  *
  */
  class InteractionHandler : public Base {
  public:
    //! Constructor
    InteractionHandler(GFlow *);

    //! @brief Add a pair of interacting particles.
    virtual void addPair(const int, const int) = 0;

    //! @brief Clear out all the data in the handler, allowing it to get new information.
    virtual void clear() = 0;
  };

}
#endif // __INTERACTION_HANDLER_HPP__GFLOW__