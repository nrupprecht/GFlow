#ifndef __MIN_INTERACTING_DISTANCE_HPP__GFLOW__
#define __MIN_INTERACTING_DISTANCE_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class MinInteractingDistance : public GraphObject {
  public:
    // Constructor
    MinInteractingDistance(GFlow*);

    virtual void post_step() override;

  private:
    //! A function we can use to compute the memory distance
    static void minimumDistanceKernel(RealType*, const RealType, const int, const int, class SimData*, const RealType*, RealType*, int);
  };

}
#endif // __MIN_INTERACTING_DISTANCE_HPP__GFLOW__