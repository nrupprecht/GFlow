#include "total-energy-data.hpp"
// Other files
#include "../../base/forcemaster.hpp"
#include "../../base/bonded.hpp"

namespace GFlowSimulation {
  // Constructor
  TotalEnergyData::TotalEnergyData(GFlow *gflow, bool ave) : GraphObject(gflow, "Energy", "time", "total energy"), useAve(ave), fraction(1.25), restrict_energy(false) {};

  void TotalEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Compute kinetic energy.
    RealType energy = 0;
    RealType **v = simData->V();
    RealType *im = simData->Im();
    int    *type = simData->Type();
    int size = Base::simData->size_owned();
    int count = 0;
    for (int n=0; n<size; ++n) {
      if (im[n]>0 && -1<type[n] && simData->isReal(n)) {
        energy += sqr(v[n], sim_dimensions)/im[n];
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
      RealType *om = simData->ScalarData(om_add), en = 0;
      RealType *sg = simData->Sg();
      for (int n=0; n<size; ++n) {
        if (im[n]>0 && -1<type[n] && simData->isReal(n)) {
          RealType II = 0.5*sqr(sg[n])/im[n];
          energy += II*sqr(om[n]);
        }
      }
      energy += 0.5*en;
    }
    
    // If this is a parallel job, sum energy from all processors.
    #if USE_MPI == 1
    MPIObject::mpi_sum(energy);
    if (useAve) MPIObject::mpi_sum(count);
    #endif

    // If we want the average
    if (useAve && count>0) energy /= static_cast<RealType>(count);

    // Store data
    RealType time = gflow->getElapsedTime();
    if (topology->getRank()==0) data.push_back(RPair(time, energy));
    // A useful check
    if(isnan(energy)) throw NanValue("Energy");
    // Set initial energy
    if (!initial_set) {
      initial_energy = energy;
      initial_set = true;
    }
    // Restrict energy
    if (restrict_energy && initial_energy>0 && energy/initial_energy>fraction) gflow->terminate();
  }

}