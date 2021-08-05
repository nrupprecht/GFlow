#include <dataobjects/multigraphobjects/groupangular.hpp>

using namespace GFlowSimulation;

GroupAngular::GroupAngular(GFlow *gflow)
    : MultiGraphObject(gflow, "GroupAngular", "time", "Mom. of Inertia", 3), Group(gflow) {
  axis_y[1] = "Ang. Momentum";
  axis_y[2] = "Torque";
};

GroupAngular::GroupAngular(GFlow *, Group &g)
    : MultiGraphObject(gflow, "GroupAngular", "time", "Mom. of Inertia", 3), Group(g) {
  axis_y[1] = "Ang. Momentum";
  axis_y[2] = "Torque";
};

void GroupAngular::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }

  // Get data.
  calculate_angular_quantities();
  // Store data.
  addEntry();
  getX() = Base::gflow->getElapsedTime();
  getY(0) = II;
  getY(1) = L;
  getY(2) = T;
}

void GroupAngular::calculate_angular_quantities() {
  if (simData == nullptr) {
    return;
  }
  // Update local ids?
  if (locals_changed) {
    update_local_ids();
  }

  // Find center of mass of the group.
  Vec com(sim_dimensions);
  findCenterOfMass(com.data);

  // Compute net angular velocity.
  Vec X(sim_dimensions), V(sim_dimensions), F(sim_dimensions);
  II = L = T = 0;

  for (int i = 0; i < Group::size(); ++i) {
    int id = at(i);
    // Get the position and force.
    X = simData->X(id);
    X -= com;
    gflow->minimumImage(X.data);

    RealType mass = 1. / simData->Im(id);
    II += mass * sqr(X); // m R^2

    V = simData->V(id); // Set Acc to be velocity.
    L += crossVec2(X, mass * V); // R x P

    F = simData->F(id);
    T += crossVec2(X, F); // R x F
  }
}

RealType GroupAngular::getII() const {
  return II;
}

RealType GroupAngular::getL() const {
  return L;
}

RealType GroupAngular::getTorque() const {
  return T;
}
