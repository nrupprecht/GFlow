#include "groupangular.hpp"

namespace GFlowSimulation {

  GroupAngular::GroupAngular(GFlow *gflow) : MultiGraphObject(gflow, "GroupAngular", "time", "Mom. of Inertia", 3) {
    axis_y[1] = "Ang. Momentum";
    axis_y[2] = "Torque";
  };

  GroupAngular::GroupAngular(GFlow*, Group& g) : MultiGraphObject(gflow, "GroupAngular", "time", "Mom. of Inertia", 3), group(g) {
    axis_y[1] = "Ang. Momentum";
    axis_y[2] = "Torque";
  };

  void GroupAngular::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get data.
    calculate_angular_quantities();
    // Store data.
    addEntry();
    getX()  = Base::gflow->getElapsedTime();
    getY(0) = II;
    getY(1) = L;
    getY(2) = T;
  }

  void GroupAngular::setGroup(Group& g) {
    group = g;
  }

  void GroupAngular::calculate_angular_quantities() {
    if (simData==nullptr) return;
    // Update local ids?
    if (locals_changed) group.update_local_ids(simData);

    cout << "Here" << endl;

    cout << simData << endl;

    // Find center of mass of the group.
    Vec com(sim_dimensions);
    group.findCenterOfMass(com.data, simData);

    cout << "Here" << endl;

    // Compute net angular velocity.
    Vec X(sim_dimensions), V(sim_dimensions), F(sim_dimensions);
    II = L = T = 0;

    cout << "Here" << endl;

    for (int i=0; i<group.size(); ++i) {

      cout << "There" << endl;
      cout << group.size() << endl;
      cout << group.at(i) << endl;

      int id = group.at(i);

      cout << id << endl;

      cout << simData->size() << endl;

      // Get the position and force.
      X = simData->X(id);
      X -= com;
      gflow->minimumImage(X.data);
      
      RealType mass = 1./simData->Im(id);
      II += mass*sqr(X); // m R^2
      
      V = simData->V(id); // Set Acc to be velocity.
      L += crossVec2(X, mass*V); // R x P

      F = simData->F(id);
      T += crossVec2(X, F); // R x F
    }
  }

  RealType GroupAngular::getII() {
    return II;
  }

  RealType GroupAngular::getL() {
    return L;
  }

  RealType GroupAngular::getTorque() {
    return T;
  }

}