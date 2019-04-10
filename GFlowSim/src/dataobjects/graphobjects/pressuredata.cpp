#include "pressuredata.hpp"
// Other files
#include "../../base/interaction.hpp"
#include "kineticenergydata.hpp"

namespace GFlowSimulation {

  PressureData::PressureData(GFlow *gflow) : GraphObject(gflow, "Pressure") {};

  void PressureData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    
    // We will look at the L1 distance between particles that potentially interact
    vector<RealType> pressures;
    // Store data
    RealType time = Base::gflow->getElapsedTime();
    RealType Ptot = calculate_pressure(gflow);
    // Store data
    data.push_back(RPair(time, Ptot));
  }

  RealType PressureData::calculate_pressure(GFlow *gflow) {
    // SimData pointer
    SimData *simData = gflow->getSimData();
    int sim_dimensions = simData->getSimDimensions();
    // Get data
    RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
    RealType KB = gflow->getKB();
    RealType T = 2.*ke/(sim_dimensions*KB);
    RealType V = gflow->getBounds().vol();
    int number = simData->number();
    RealType factor = 1./(sim_dimensions*V);
    RealType Ptot = number*KB*T/V;
    const vector<Interaction*>& interactions = gflow->getInteractions();
    // Get the virials from all the interactions
    for (const auto it : interactions) {
      //! P = N k T/V + 1/(sim_dimensions*V) \sum_i (r_i \dot F_i)
      // virial = \sum_i (r_i \dot F_i)
      RealType virial = it->getVirial();
      // Update Ptot.
      Ptot += factor*virial;
    }

    // Normalize
    Ptot /= (number*KB*T/V);

    // Return
    return Ptot;
  }


}