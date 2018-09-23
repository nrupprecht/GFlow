#ifndef __MIN_INTERACTING_DISTANCE_HPP__GFLOW__
#define __MIN_INTERACTING_DISTANCE_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class MinInteractingDistance : public DataObject {
  public:
    // Constructor
    MinInteractingDistance(GFlow*);

    virtual void post_step() final;

    virtual bool writeToFile(string, bool=true) final;

  private:
    vector<RPair> data;

    //! A function we can use to compute the memory distance
    static void minimumDistanceKernel(RealType*, const RealType, const int, const int, class SimData*, const RealType*, RealType*);
  };

}
#endif // __MIN_INTERACTING_DISTANCE_HPP__GFLOW__