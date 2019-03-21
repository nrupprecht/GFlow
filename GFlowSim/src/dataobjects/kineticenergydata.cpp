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
        ke += sqr(v[n], sim_dimensions)/im[n];
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

}