#include "integrator.hpp"

namespace GFlowSimulation {

  Integrator::Integrator(GFlow *gflow) : Base(gflow), dt(0.0001) {};

  RealType Integrator::getTimeStep() {
    return dt;
  }

}