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
    RealType *om = simData->ScalarData(om_add);
    RealType *im = simData->Im();
    RealType *sg = simData->Sg();
    int    *type = simData->Type();
    int size = Base::simData->size();
    int count = 0;

    // Check if there is rotational motion.
    if (om==nullptr) return;

    for (int n=0; n<size; ++n) {
      if (im[n]>0 && -1<type[n] && simData->isReal(n)) {
        RealType II = 0.5*sqr(sg[n])/im[n];
        energy += II*sqr(om[n]);
        ++count;
      }
    }
    energy *= 0.5;

    // Store data
    RealType time = gflow->getElapsedTime();
    data.push_back(RPair(time, energy));
    // A useful check
    if(isnan(energy)) throw NanValue("Energy");
  }

}