#ifndef __TIME_STEP_DATA_HPP__GFLOW__
#define __TIME_STEP_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class TimeStepData : public GraphObject {
  public:
    // Constructor
    TimeStepData(GFlow*);

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;
  };

}
#endif // __TIME_STEP_DATA_HPP__GFLOW__