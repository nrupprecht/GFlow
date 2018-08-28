#ifndef __MEMORY_DISTANCE_HPP__GFLOW__
#define __MEMORY_DISTANCE_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class MemoryDistance : public DataObject {
  public:
    //! Constructor
    MemoryDistance(GFlow*);

    //! Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //! Returns true for success.
    virtual bool writeToFile(string, bool=true);

    //! @brief This executes the statistic, finding the average distance in memory between potentially
    //! interacting particles.
    RealType getAverageDistance();

  private:
    //! The vector of data.
    vector<RPair> data;

    //! A function we can use to compute the memory distance
    static void memoryDistanceKernel(RealType*, const RealType, const int, const int, class SimData*, const RealType*, RealType*);
  };

}
#endif // __MEMORY_DISTANCE_HPP__GFLOW__