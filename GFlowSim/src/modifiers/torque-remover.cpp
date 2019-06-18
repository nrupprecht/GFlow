#include "torque-remover.hpp"
// Other files
#include "../dataobjects/graphobjects/grouptorque.hpp"

namespace GFlowSimulation {

  TorqueRemover::TorqueRemover(GFlow *gflow) : Modifier(gflow) {};

  TorqueRemover::TorqueRemover(GFlow *gflow, Group& g) : Modifier(gflow), group(g) {};

  void TorqueRemover::post_forces() {
    Modifier::post_forces();

    // We need there to be a valid group.
    if (group.empty()) return;

    // This only works (for now?) in two dimensions.
    if (sim_dimensions!=2) return;

    // Update local ids?
    if (simData->getNeedsRemake()) group.update_local_ids(simData);

    // Find center of mass of the group.
    Vec com(sim_dimensions);
    group.findCenterOfMass(com.data, simData);

    // Compute total torque. We must use torque = dII/dt * omega + II * alpha 
    // since II is non-constant.
    Vec X(sim_dimensions), F(sim_dimensions);
    RealType omega = 0, torque = 0, II = 0;
    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);

      // Get the position and force.
      X = simData->X(id);
      X -= com;
      gflow->minimumImage(X.data);
      F = simData->F(id);
      
      RealType mass = 1./simData->Im(id);
      II += mass*sqr(X);
      // Update torque. Tq = R x F
      torque += crossVec2(X.data, F.data); 
    }

    // Angular acceleration
    RealType alpha = II>0 ? torque/II : 0;

    cout << torque << ", " << II << endl;

    Vec Fnet(sim_dimensions);

    // Apply torque to counteract the angular acceleration.
    for (int i=0; i<group.size(); ++i) {
      int id = group.at(i);
      
      // Calculate requisite force.
      X = simData->X(id);
      X -= com;
      gflow->minimumImage(X.data);
      //gflow->getDisplacement(simData->X(id), com.data, X.data);

      RealType mass = 1./simData->Im(id);
      F[0] =   mass * alpha * X[1];
      F[1] =  -mass * alpha * X[0];

      //cout << alpha << " " << F << " " << " {" << simData->F(id, 0) <<  ", "  << simData->F(id, 0) << "}\n";

      simData->F(id, 0) += F[0];
      simData->F(id, 1) += F[1];

      Fnet += F;
    }
    
    //cout << "Torque: " << GroupTorque::calculate_torque(simData, group) << ", Fnet = " << Fnet <<  endl;
  }

  void TorqueRemover::setGroup(Group& g) {
    group = g;
  }

}