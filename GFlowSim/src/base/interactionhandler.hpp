#ifndef __INTERACTION_HANDLER_HPP__GFLOW__
#define __INTERACTION_HANDLER_HPP__GFLOW__

#include "../gflow.hpp"
#include "../utility/vectormath.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  //! @brief This defines what type of force kernel function the interaction handler expects to be passed in to it.
  // using Kernel = auto (*) (simd_float*, simd_float*, const simd_float, const simd_float, const simd_float*, const RealType*, RealType*) -> void;

  //template<typename T>
  //using Kernel = auto (*) (T*, T*, const T, const T, const T*, const RealType*, RealType*) -> void;

  template<typename float_type>
  using Kernel = auto (*) (float_type*, const float_type*, const float_type, const float_type, const float_type*, const RealType*, RealType *) -> void;

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
    InteractionHandler(GFlow *gflow) : Base(gflow) {};

    //! @brief Add a pair of interacting particles.
    virtual void addPair(const int, const int) = 0;

    //! @brief Signals that the pair additions are done.
    virtual void close() = 0;

    //! @brief Clear out all the data in the handler, allowing it to get new information.
    virtual void clear() = 0;

    //! @brief Returns a characteriztic size indicating the number of interactions the handler will test for. 
    //! 
    //! This may be the size of the array the particles are stored in, or the number of interactions the
    //! handler will try. If unknown, the function returns -1. This could occur if e.g. we are travering 
    //! through linked cells, and don't have an accurate guess of how many interactions there will be.
    virtual int size() const = 0;

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
    virtual void executeKernel(Kernel<simd_float>, Kernel<float>, const RealType*, RealType*, const vector<int>&) const = 0;
  };

}
#endif // __INTERACTION_HANDLER_HPP__GFLOW__