#include "rotational-energy-data.hpp"
// Other files
#include "../../base/forcemaster.hpp"
#include "../../base/bonded.hpp"

namespace GFlowSimulation {
  // Constructor
  RotationalEnergyData::RotationalEnergyData(GFlow *gflow, bool ave) : GraphObject(gflow, "RotationalEnergy", "time", "rotational energy"), useAve(ave) {};

  void RotationalEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Compute kinetic energy.
    RealType energy = 0;
    int om_add = simData->getScalarData("Om");
    auto om = simData->ScalarData(om_add);
    auto im = simData->Im();
    auto sg = simData->Sg();
    auto type = simData->Type();
    int count = 0;

    // Check if there is rotational motion.
    if (om.isnull()) return;

    for (int n=0; n<simData->size_owned(); ++n) {
      if (im[n]>0 && -1<type[n]) {
        RealType II = 0.5*sqr(sg[n])/im[n];
        energy += II*sqr(om[n]);
        ++count;
      }
    }
    energy *= 0.5;

    // Store data. These functions work correctly with multiprocessor runs.
    if (useAve) gatherAverageData(gflow->getElapsedTime(), energy, count);
    else gatherData(gflow->getElapsedTime(), energy);
  }

}