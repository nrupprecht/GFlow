#include "total-energy-data.hpp"
// Other files
#include "../../base/forcemaster.hpp"
#include "../../base/bonded.hpp"

namespace GFlowSimulation {
  // Constructor
  TotalEnergyData::TotalEnergyData(GFlow *gflow, bool ave) : GraphObject(gflow, "Energy", "time", "total energy"), useAve(ave) {};

  void TotalEnergyData::pre_integrate() {
    forceMaster->setCalculatePotential(true);
  }

  void TotalEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Compute kinetic energy.
    RealType energy = 0;
    auto v = simData->V();
    auto im = simData->Im();
    auto type = simData->Type();
    int size = simData->size_owned();
    int count = 0;
    for (int n=0; n<size; ++n) {
      if (im[n]>0 && -1<type[n] && simData->isReal(n)) {
        energy += sqr(v(n), sim_dimensions)/im[n];
        ++count;
      }
    }
    energy *= 0.5;
    // Add potential energy. \todo This may not be accurate for multicore runs, some potential may be double counted.
    energy += forceMaster->getTotalPotentialEnergy();
    // Get boundary energy from gflow.
    energy += gflow->getBoundaryEnergy();
    // Add potential energy from bonded interactions
    auto bonded_interactions = gflow->getBondedInteractions();
    for (auto it : bonded_interactions) energy += it->getPotential();
    // Look for rotational energy
    int om_add = simData->getScalarData("Om");
    if (om_add>=0) {
      RealType en = 0;
      auto om = simData->ScalarData(om_add);
      auto sg = simData->Sg();
      for (int n=0; n<size; ++n) {
        if (im[n]>0 && -1<type[n] && simData->isReal(n)) {
          RealType II = 0.5*sqr(sg[n])/im[n];
          energy += II*sqr(om[n]);
        }
      }
      energy += 0.5*en;
    }

    // Store data. These functions work correctly with multiprocessor runs.
    if (useAve) gatherAverageData(gflow->getElapsedTime(), energy, count);
    else gatherData(gflow->getElapsedTime(), energy);

  }

}