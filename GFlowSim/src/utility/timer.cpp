#include "timer.hpp"

namespace GFlowSimulation {

  void Timer::start() {
    _s = current_time();
  }
  
  void Timer::stop() {
    auto _e = current_time();
    _t += time_span(_e, _s);
  }

  void Timer::clear() {
    _t = 0.;
  }

  double Timer::time() {
    return _t;
  }


}