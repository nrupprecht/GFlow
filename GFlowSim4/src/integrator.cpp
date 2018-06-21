#include "integrator.hpp"

namespace GFlowSimulation {

  Integrator::Integrator(GFlow *gflow) : Base(gflow), dt(DEFAULT_TIME_STEP) {};

  RealType Integrator::getTimeStep() {
    return dt;
  }

}