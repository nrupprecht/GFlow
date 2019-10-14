#include "kineticenergydata.hpp"

namespace GFlowSimulation {
  // Constructor
  KineticEnergyData::KineticEnergyData(GFlow *gflow, bool ave) : GraphObject(gflow, "KE", "time", "kinetic energy"), useAve(ave) {};

  void KineticEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // Get data.
    RealType ke = calculate_kinetic(simData, false);
    // Store data. These functions work correctly with multiprocessor runs.
    if (useAve) gatherAverageData(gflow->getElapsedTime(), ke, simData->size());
    else gatherData(gflow->getElapsedTime(), ke);
  }

  RealType KineticEnergyData::calculate_kinetic(shared_ptr<SimData> simData, bool average) {
    RealType ke = 0;
    auto v = simData->V();
    auto im = simData->Im();
    int size = simData->size();
    int sim_dimensions = simData->getSimDimensions();
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && simData->isReal(n)) {
        RealType m = 1./im[n];
        ke += m*sqr(v[n], sim_dimensions);
        ++count;
      }
    ke *= 0.5;

    // Return the total kinetic energy
    return average ? ke/static_cast<RealType>(count) : ke;
  }

  RealType KineticEnergyData::calculate_temperature(shared_ptr<SimData> simData) {
    // Get the kinetic energy.
    RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
    // Calculate temperature.
    return 2./simData->getSimDimensions() * ke / simData->getGFlow()->getKB();
  }

}