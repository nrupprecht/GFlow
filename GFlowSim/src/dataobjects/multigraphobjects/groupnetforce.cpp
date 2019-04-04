#include "groupnetforce.hpp"

namespace GFlowSimulation {

  GroupNetForce::GroupNetForce(GFlow *gflow) 
    : MultiGraphObject(gflow, "GroupForce", "time", "force", gflow->getSimDimensions()) {};

  void GroupNetForce::post_step() {
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

    // Set the forces
    RealType *F = new RealType[sim_dimensions];
    group.findNetForce(F, simData);
    // Set each dimension of force.
    for (int d=0; d<sim_dimensions; ++d) getY(d) = F[d];

    // Clean up
    delete [] F;
  }

  void GroupNetForce::setGroup(const Group& g) {
    group = g;
  }

}