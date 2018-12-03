#ifndef __TIMER_HPP__GFLOW__
#define __TIMER_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  class Timer {
  public:
    void start();
    void stop();
    void clear();
    double time();
  private:
    high_resolution_clock::time_point _s;
    double _t = 0;
  };

}
#endif // __TIMER_HPP__GFLOW__