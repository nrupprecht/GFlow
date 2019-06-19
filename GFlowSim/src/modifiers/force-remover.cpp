#include "force-remover.hpp"

namespace GFlowSimulation {

  ForceRemover::ForceRemover(GFlow *gflow) : Modifier(gflow), projection(gflow->getSimDimensions()) {
    // Default value is y-hat;
    projection[1] = 1.;
  };

  ForceRemover::ForceRemover(GFlow *gflow, RealType *v) : Modifier(gflow), projection(gflow->getSimDimensions()) {
    copyVec(v, projection);
  }

  ForceRemover::ForceRemover(GFlow *gflow, Vec& v) : Modifier(gflow), projection(v) {};

  ForceRemover::ForceRemover(GFlow *gflow, Group& group) : Modifier(gflow), projection(gflow->getSimDimensions()) {
    // Default value is y-hat;
    projection[1] = 1.;
    // Set the group.
    set(group);
  }

  void ForceRemover::pre_integrate() {
    Modifier::pre_integrate();

    // Update local ids
    update_local_ids(simData);

    // Find the total mass of the group.
    mass = findTotalMass(simData);
    // Make sure projection vector is normalized.
    projection.normalize();

    // Remove any velocity in the [projection] direction
    Vec V(sim_dimensions);
    findCOMVelocity(V.data, simData);
    V = -(V*projection)*projection;
    addVelocity(V.data, simData);
    // Remove any initial force (which should be zero anyways).
    remove_net_force();
  }

  void ForceRemover::post_forces() {
    Modifier::post_forces();
    // Remove the net force.
    remove_net_force();
  }

  inline void ForceRemover::remove_net_force() {
    if (mass<=0 || empty()) return;

    // Update local ids
    if (simData->getNeedsRemake()) update_local_ids(simData);

    // Find the net acceleration of the group of particles.
    Vec F(sim_dimensions);
    findNetForce(F.data, simData);
    RealType acceleration = (projection*F)/mass;
    // Turn F into an acceleration vector.
    F = -acceleration*projection;
    // Add an acceleration to cancel out the net acceleration of the group.
    addAcceleration(F.data, simData);
  }

}