#include "grouptorque.hpp"

namespace GFlowSimulation {

  GroupTorque::GroupTorque(GFlow *gflow) : GraphObject(gflow, "Torque", "time", "Torque"), Group(gflow) {};

  GroupTorque::GroupTorque(GFlow*, Group& g) : GraphObject(gflow, "Torque", "time", "Torque"), Group(g) {};

  void GroupTorque::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // Update local ids?
    if (locals_changed) update_local_ids();
    locals_changed = false;
    // Store data. These functions work correctly with multiprocessor runs.
    gatherData(gflow->getElapsedTime(), calculate_torque(simData, *static_cast<Group*>(this)));
  }

  RealType GroupTorque::calculate_torque(shared_ptr<SimData> simData, const Group& group) {
    // Get gflow and the number of dimensions.
    int sim_dimensions = simData->getSimDimensions();
    GFlow *gflow = simData->getGFlow();

    // Find center of mass of the group.
    Vec com(sim_dimensions);
    group.findCenterOfMass(com.data);

    // Compute total torque.
    Vec X(sim_dimensions), Acc(sim_dimensions);
    RealType omega = 0, torque = 0, II = 0, dII = 0;
    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);

      // Get the position and force.
      X = simData->X(id);
      Acc = simData->F(id);
      X -= com;
      gflow->minimumImage(X.data);
      torque += crossVec2(X.data, Acc.data);
    }
    return torque;
  }

}