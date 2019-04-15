#include "timedobject.hpp"

namespace GFlowSimulation {

  void TimedObject::start_timer() {
    timer.start();
  }

  void TimedObject::stop_timer() {
    timer.stop();
  }

  void TimedObject::clear_timer() {
    timer.clear();
  }

  double TimedObject::get_time() {
    return timer.time();
  }

}