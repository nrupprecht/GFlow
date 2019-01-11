#include "langevintypeintegrator.hpp"

namespace GFlowSimulation {

  LangevinTypeIntegrator::LangevinTypeIntegrator(GFlow *gflow, RealType T, RealType vis) : Integrator(gflow), viscosity(vis), 
  temperature(T), lastUpdate(0), updateDelay(0.01)
  {
    drift1 = temperature/(6.*viscosity*PI);
  }

  void LangevinTypeIntegrator::setViscosity(RealType eta) {
    viscosity = eta;
    drift1 = temperature/(6.*viscosity*PI);
  }

  void LangevinTypeIntegrator::setTemperature(RealType T) {
    temperature = T;
    drift1 = temperature/(6.*viscosity*PI);
  }

}
