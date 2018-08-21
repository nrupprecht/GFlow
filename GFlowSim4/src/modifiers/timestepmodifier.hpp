#ifndef __TIME_STEP_MODIFIER_HPP__GFLOW__
#define __TIME_STEP_MODIFIER_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class TimestepModifier : public Modifier {
  public:
    // Constructor
    TimestepModifier(GFlow*);

    virtual void pre_integrate() final;
    virtual void post_forces() final;

  private:
    RealType minSigma, vtollerance, atollerance;
    RealType lastCheck, delay;
    RealType max_dt, min_dt;
  };

}
#endif // __TIME_STEP_MODIFIER_HPP__GFLOW__