#include "bonded-energy-data.hpp"
// Other files
#include "../../base/forcemaster.hpp"
#include "../../base/bonded.hpp"

namespace GFlowSimulation {
  // Constructor
  BondedEnergyData::BondedEnergyData(GFlow *gflow, bool ave) : GraphObject(gflow, "BondEnergy", "time", "bond energy"), useAve(ave) {};

  void BondedEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType energy = 0;
    // Add potential energy from bonded interactions
    auto bonded_interactions = gflow->getBondedInteractions();
    for (auto it : bonded_interactions) energy += it->getPotential();
    
    // Store data. These functions work correctly with multiprocessor runs.
    if (useAve) gatherAverageData(gflow->getElapsedTime(), energy, simData->size_owned());
    else gatherData(gflow->getElapsedTime(), energy);
  }

}