#include "groupnetforce.hpp"

namespace GFlowSimulation {

  GroupNetForce::GroupNetForce(GFlow *gflow) : MultiGraphObject(gflow, "GroupForce", "time", "force", gflow->getSimDimensions()), Group(gflow) {};

  void GroupNetForce::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // No particles to keep track of.
    if (size()==0) return;

    // Make sure local ids are up to date.
    if (locals_changed) update_local_ids();
    locals_changed = false;

    // Add a new entry to modify
    addEntry();

    // Set the time
    getX() = Base::gflow->getElapsedTime();
    // Set the forces
    Vec F(sim_dimensions);
    findNetForce(F.data);
    // Set each dimension of force.
    for (int d=0; d<sim_dimensions; ++d) getY(d) = F[d];
  }

}