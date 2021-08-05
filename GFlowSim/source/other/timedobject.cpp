#include "other/timedobject.hpp"

using namespace GFlowSimulation;

// Initialize timing on to be true.
bool TimedObject::timing_on = true;

void TimedObject::start_timer() {
  if (TimedObject::timing_on) {
    timer.start();
  }
}

void TimedObject::stop_timer() {
  if (TimedObject::timing_on) {
    timer.stop();
  }
}

void TimedObject::clear_timer() {
  timer.clear();
}

double TimedObject::get_time() {
  return timer.time();
}

bool TimedObject::getTimingOn() {
  return TimedObject::timing_on;
}

void TimedObject::setTimingOn(bool t) {
  TimedObject::timing_on = t;
}
