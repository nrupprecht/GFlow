#include <dataobjects/multigraphobjects/groupnetforce.hpp>

using namespace GFlowSimulation;

GroupNetForce::GroupNetForce(GFlow *gflow)
    : MultiGraphObject(gflow, "GroupForce", "time", "force", gflow->getSimDimensions()), Group(gflow) {
  for (int i = 0; i < sim_dimensions; ++i) {
    axis_y[i] = "Force - F[" + toStr(i) + "]";
  }
};

void GroupNetForce::post_step() {
  // Only record if enough time has gone by and there are particles to keep track of.
  if (!DataObject::_check()) {
    return;
  }

  // Make sure local ids are up to date.
  if (locals_changed) {
    update_local_ids();
    locals_changed = false;
  }

  // Gather and store data on processor 0.
  Vec F(sim_dimensions);
  findNetForce(F.data);
  gatherData(gflow->getElapsedTime(), F);
}
