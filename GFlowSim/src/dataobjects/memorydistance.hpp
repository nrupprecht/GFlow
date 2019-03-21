#ifndef __MEMORY_DISTANCE_HPP__GFLOW__
#define __MEMORY_DISTANCE_HPP__GFLOW__

#include "graphobject.hpp"

namespace GFlowSimulation {

  class MemoryDistance : public GraphObject {
  public:
    //! Constructor
    MemoryDistance(GFlow*);

    //! Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! @brief This executes the statistic, finding the average distance in memory between potentially
    //! interacting particles.
    RealType getAverageDistance();

  private:
    //! A function we can use to compute the memory distance
    static void memoryDistanceKernel(RealType*, const RealType, const int, const int, class SimData*, const RealType*, RealType*);
  };

}
#endif // __MEMORY_DISTANCE_HPP__GFLOW__