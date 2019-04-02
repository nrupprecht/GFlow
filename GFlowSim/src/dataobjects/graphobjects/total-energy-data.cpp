#include "total-energy-data.hpp"
// Other files
#include "../../base/forcemaster.hpp"

namespace GFlowSimulation {
  // Constructor
  TotalEnergyData::TotalEnergyData(GFlow *gflow, bool ave) : GraphObject(gflow, "Energy", "time", "total energy"), useAve(ave) {};

  void TotalEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType energy = 0;
    RealType **v = simData->V();
    RealType *im = simData->Im();
    int    *type = simData->Type();
    int size = Base::simData->size();
    int count = 0;
    for (int n=0; n<size; ++n) {
      if (im[n]>0 && -1<type[n]) {
        energy += sqr(v[n], sim_dimensions)/im[n];
        ++count;
      }
    }
    energy *= 0.5;
    // Add potential energy
    energy += forceMaster->getTotalPotentialEnergy();
    // If we want the average
    if (useAve && count>0) energy /= count;
    // Store data
    RealType time = gflow->getElapsedTime();
    data.push_back(RPair(time, energy));
    // A useful check
    if(isnan(energy)) throw NanValue("Energy");
  }

}