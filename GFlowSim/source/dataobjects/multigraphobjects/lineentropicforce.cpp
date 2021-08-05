#include <dataobjects/multigraphobjects/lineentropicforce.hpp>
// Other
#include <utility/vec.hpp>
#include <dataobjects/graphobjects/kineticenergydata.hpp>

using namespace GFlowSimulation;

LineEntropicForce::LineEntropicForce(GFlow *gflow)
    : MultiGraphObject(gflow, "LineEntropicForce", "time", "F/(rho kB T L)", gflow->getSimDimensions()),
      group(gflow) {};

void LineEntropicForce::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  // Make sure local ids are up to date.
  if (locals_changed) {
    group.update_local_ids();
  }
  locals_changed = false;

  // Add a new entry to modify
  addEntry();

  // Store data
  RealType time = Base::gflow->getElapsedTime();
  // Set the time
  getX() = time;

  // Compute the temperature of the system.
  int n_solvent = simData->number() - 2 * group.size();
  RealType KB = gflow->getKB();
  // Get KE per (non infinitely massive) particle.
  RealType ke = KineticEnergyData::calculate_kinetic(simData, true);
  RealType T = 2. / static_cast<RealType>(sim_dimensions) * ke / KB;
  // Compute the number density. Do not include your own particles. Assumes the other line has roughly the same number of particles.
  RealType rho = n_solvent / gflow->getBounds().vol();
  // Compute the normalization
  RealType norm = 1. / (rho * KB * T * length);

  // Set the forces.
  Vec F(sim_dimensions);
  group.findNetForce(F.data);
  // Data reduction.
  MPIObject::mpi_sum0(F.data, sim_dimensions);

  // Set each dimension of force.
  for (int d = 0; d < sim_dimensions; ++d) {
    getY(d) = norm * F[d];
  }
}

void LineEntropicForce::setLength(RealType l) {
  length = l;
}

void LineEntropicForce::setGroup(const Group &g, bool setLength) {
  group = g;
  if (setLength) {
    // Set length
    int first = simData->getLocalID(group.g_at(0));
    int last = simData->getLocalID(group.g_at(group.size() - 1));
    length = simData->X(last, 1) - simData->X(first, 1);
  }
}
