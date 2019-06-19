#include "torque-remover.hpp"
// Other files
#include "../dataobjects/graphobjects/grouptorque.hpp"
#include "../dataobjects/multigraphobjects/groupangular.hpp"

namespace GFlowSimulation {

  TorqueRemover::TorqueRemover(GFlow *gflow) : Modifier(gflow) {};

  TorqueRemover::TorqueRemover(GFlow *gflow, Group& g) : Modifier(gflow), group(g) {};

  void TorqueRemover::pre_integrate() {
    group.update_local_ids(simData);

    // Find center of mass of the group.
    Vec com(sim_dimensions), X(sim_dimensions), V(sim_dimensions);
    group.findCenterOfMass(com.data, simData);

    // Compute total antular velocity.
    RealType L = 0, II = 0;
    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);
      // Get the position and force.
      gflow->getDisplacement(simData->X(id), com.data, X.data);
      // Compute
      RealType mass = 1./simData->Im(id);
      II += mass*sqr(X);
      V = simData->V(id);
      L += crossVec2(X, mass*V);
    }
    // Calculate angular velocity
    RealType omega = L / II;

    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);
      // Remove angular velocity.
      gflow->getDisplacement(simData->X(id), com.data, X.data);
      V[0] =   omega * X[1];
      V[1] =  -omega * X[0];
      simData->V(id, 0) += V[0];
      simData->V(id, 1) += V[1];
    }
  }

  void TorqueRemover::post_forces() {
    Modifier::post_forces();

    // We need there to be a valid group.
    if (group.empty() || simData==nullptr) return;

    // This only works (for now?) in two dimensions.
    if (sim_dimensions!=2) return;

    // Update local ids?
    if (simData->getNeedsRemake()) group.update_local_ids(simData);

    // Find center of mass of the group.
    Vec com(sim_dimensions);
    group.findCenterOfMass(com.data, simData);

    // Compute total torque. We must use torque = dII/dt * omega + II * alpha 
    // since II is non-constant.
    Vec X(sim_dimensions), V(sim_dimensions), F(sim_dimensions);
    RealType L = 0, torque = 0, II = 0;
    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);
      //X = simData->X(id);
      //X -= com;
      //gflow->minimumImage(X.data);
      gflow->getDisplacement(simData->X(id), com.data, X.data);
      
      RealType mass = 1./simData->Im(id);
      II += mass*sqr(X);
      // Update torque. Tq = R x F
      torque += crossVec2(X.data, simData->F(id)); 
      L += crossVec2(X, mass*V);
    }

    // Angular acceleration
    RealType alpha = II>0 ? torque/II : 0;

    // Apply torque to counteract the angular acceleration.
    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);
      gflow->getDisplacement(simData->X(id), com.data, X.data);
      RealType mass = 1./simData->Im(id);
      simData->F(id, 0) += mass * alpha * X[1];
      simData->F(id, 1) -= mass * alpha * X[0];
    }
  }

  void TorqueRemover::setGroup(Group& g) {
    group = g;
  }

}