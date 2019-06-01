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
    // If we want the average
    if (useAve && Base::simData->size()>0) energy /= Base::simData->size();;
    // Store data
    RealType time = gflow->getElapsedTime();
    data.push_back(RPair(time, energy));
    // A useful check
    if(isnan(energy)) throw NanValue("Bonded Energy");
  }

}