#include "average-omega-data.hpp"

namespace GFlowSimulation {

  // Constructor
  AverageOmegaData::AverageOmegaData(GFlow *gflow) : GraphObject(gflow, "AverageOmega", "time", "omega") {};

  void AverageOmegaData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Compute kinetic energy.
    int om_add = simData->getScalarData("Om");
    auto om = simData->ScalarData(om_add);
    auto im = simData->Im();
    auto type = simData->Type();
    int size = Base::simData->size();

    // Check if there is rotational motion.
    if (om.isnull()) return;

    // Average the angular velocity, weighing by particle mass.
    RealType mass = 0, omega = 0;
    for (int n=0; n<size; ++n) {
      if (im[n]>0 && -1<type[n] && simData->isReal(n)) {
        RealType m = 1./im[n];
        omega += om[n]*m;
        mass += m;
      }
    }

    // Store data. These functions work correctly with multiprocessor runs.
    gatherAverageData(gflow->getElapsedTime(), omega, mass);
  }

}