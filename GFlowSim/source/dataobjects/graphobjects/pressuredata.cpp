#include <dataobjects/graphobjects/pressuredata.hpp>
// Other files
#include <base/interaction.hpp>
#include <dataobjects/graphobjects/kineticenergydata.hpp>

using namespace GFlowSimulation;

PressureData::PressureData(GFlow *gflow)
    : GraphObject(gflow, "Pressure") {};

void PressureData::pre_integrate() {
  const auto &interactions = gflow->getInteractions();
  // Get the virials from all the interactions
  for (const auto it : interactions) {
    it->setDoVirial(true);
  }
}

void PressureData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  // Calculate pressure.
  RealType total_pressure = calculate_pressure(gflow);

  // Store data. These functions work correctly with multiprocessor runs.
  gatherData(gflow->getElapsedTime(), total_pressure);
}

RealType PressureData::calculate_pressure(GFlow *gflow) {
  // SimData pointer
  auto simData = gflow->getSimData();
  int sim_dimensions = simData->getSimDimensions();
  // Get data
  RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
  RealType KB = gflow->getKB();
  RealType T = 2. * ke / (sim_dimensions * KB);
  RealType V = gflow->getBounds().vol();
  int number = simData->number();
  RealType factor = 1. / (sim_dimensions * V);
  RealType Ptot = number * KB * T / V;
  const auto &interactions = gflow->getInteractions();
  // Get the virials from all the interactions
  for (const auto it : interactions) {
    RealType virial = it->getVirial();
    Ptot += factor * virial;
  }

  // Normalize
  Ptot /= (number * KB * T / V);

  // Return
  return Ptot;
}
