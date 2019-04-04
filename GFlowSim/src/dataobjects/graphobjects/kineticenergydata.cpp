#include "kineticenergydata.hpp"

namespace GFlowSimulation {
  // Constructor
  KineticEnergyData::KineticEnergyData(GFlow *gflow, bool ave) : GraphObject(gflow, "KE", "time", "kinetic energy"), useAve(ave) {};

  void KineticEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType ke = 0;
    RealType **v = Base::simData->V();
    RealType *im = Base::simData->Im();
    int size = Base::simData->size();
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0) {
        RealType m = 1./im[n];
        ke += m*sqr(v[n], sim_dimensions);
        ++count;
      }
    ke *= 0.5;
    // If we want the average
    if (useAve && count>0) ke /= count;
    // Store data
    RealType time = Base::gflow->getElapsedTime();
    data.push_back(RPair(time, ke));

    // A useful check
    if(isnan(ke)) throw NanValue("KE");
  }

  RealType KineticEnergyData::calculate_kinetic(SimData *simData, bool average) {
    RealType ke = 0;
    RealType **v = simData->V();
    RealType *im = simData->Im();
    int size = simData->size();
    int sim_dimensions = simData->getSimDimensions();
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0) {
        RealType m = 1./im[n];
        ke += m*sqr(v[n], sim_dimensions);
        ++count;
      }
    ke *= 0.5;
    // Return the total kinetic energy
    return average ? ke/static_cast<RealType>(count) : ke;
  }

}