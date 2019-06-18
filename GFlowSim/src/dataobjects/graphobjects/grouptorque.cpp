#include "grouptorque.hpp"

namespace GFlowSimulation {

  GroupTorque::GroupTorque(GFlow *gflow) : GraphObject(gflow, "Torque", "time", "Torque") {};

  GroupTorque::GroupTorque(GFlow*, Group& g) : GraphObject(gflow, "Torque", "time", "Torque"), group(g) {};

  void GroupTorque::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get data.
    RealType time = Base::gflow->getElapsedTime();
    RealType torque = calculate_torque(simData, group);
    // Store data.
    data.push_back(RPair(time, torque));
  }

  void GroupTorque::setGroup(Group& g) {
    group = g;
  }

  RealType GroupTorque::calculate_torque(SimData *simData, const Group& group) {
    // Update local ids?
    if (simData->getNeedsRemake()) group.update_local_ids(simData);

    // Get gflow and the number of dimensions.
    int sim_dimensions = simData->getSimDimensions();
    GFlow *gflow = simData->getGFlow();

    // Find center of mass of the group.
    Vec com(sim_dimensions);
    group.findCenterOfMass(com.data, simData);

    // Compute total torque. We must use torque = dII/dt * omega + II * alpha 
    // since II is non-constant.
    Vec X(sim_dimensions), Acc(sim_dimensions);
    RealType omega = 0, torque = 0, II = 0, dII = 0;
    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);

      // Get the position and force.
      X = simData->X(id);
      Acc = simData->F(id);
      X -= com;
      gflow->minimumImage(X.data);
      //RealType mass = 1./simData->Im(id);
      //II += mass*sqr(X);
      // Update omega, torque.
      torque += crossVec2(X.data, Acc.data);
      
      //Acc = simData->V(id); // Set Acc to be velocity.
      //RealType rsqr = sqr(X);
      //omega  += rsqr>0 ? crossVec2(X, Acc) / rsqr : 0;
      //dII += mass*(X*Acc);
    }
    //dII *= 2;
    //torque += dII*omega;

    return torque;
  }

}