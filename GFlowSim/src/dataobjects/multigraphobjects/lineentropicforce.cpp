#include "lineentropicforce.hpp"
// Other
#include "../../utility/vec.hpp"
#include "../graphobjects/kineticenergydata.hpp"

namespace GFlowSimulation {

  LineEntropicForce::LineEntropicForce(GFlow *gflow) 
    : MultiGraphObject(gflow, "LineEntropicForce", "time", "F/(rho kB T L)", gflow->getSimDimensions()) {};

  void LineEntropicForce::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // No particles to keep track of.
    if (group.size()==0) return;

    // Make sure local ids are up to date.
    if (locals_changed) group.update_local_ids(simData);
    locals_changed = false;

    // Add a new entry to modify
    addEntry();

    // Store data
    RealType time = Base::gflow->getElapsedTime();
    // Set the time
    getX() = time;

    // Compute the temperature of the system.
    RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
    RealType KB = gflow->getKB();
    RealType T = 2./static_cast<RealType>(sim_dimensions) * ke / KB;
    // Compute the number density. Do not include your own particles. Assumes the other line has roughly the same number of particles.
    RealType rho = (simData->number() - 2*group.size()) / gflow->getBounds().vol();
    // Compute the normalization
    RealType norm = 1./(rho * KB * T * length);

    // Set the forces
    Vec F(sim_dimensions);
    group.findNetForce(F.data, simData);
    // Set each dimension of force.
    for (int d=0; d<sim_dimensions; ++d) getY(d) = norm*F[d]; 
  }

  void LineEntropicForce::setLength(RealType l) {
    length = l;
  }

  void LineEntropicForce::setGroup(const Group& g) {
    group = g;
  }

}