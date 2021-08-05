#include "utility/timer.hpp"

using namespace GFlowSimulation;

void Timer::start() {
  _s = current_time();
  _running = true;
}

void Timer::stop() {
  if (_running) {
    auto _e = current_time();
    _last_t = time_span(_e, _s);
    _t += _last_t;
    _running = false;
  }
}

void Timer::clear() {
  _t = 0.;
}

double Timer::time() {
  return _t;
}

double Timer::last() {
  return _last_t;
}

double Timer::current() {
  if (_running) {
    auto _e = current_time();
    return _t + time_span(_e, _s);
  }
  else {
    return _t;
  }
}
