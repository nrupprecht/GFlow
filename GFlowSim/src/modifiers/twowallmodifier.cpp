#include "twowallmodifier.hpp"
// Other files
#include "../utility/vec.hpp"
#include "../dataobjects/graphobjects/twowallbinforce.hpp"

namespace GFlowSimulation {

  TwoWallModifier::TwoWallModifier(GFlow *gflow) : Modifier(gflow), wallA(WallSlideBody(gflow)), wallB(WallSlideBody(gflow)), 
    max_distance(0.25), acceleration(10.), dissipation(0.) {
    // Add the bodies to gflow
    gflow->addBody(&wallA);
    gflow->addBody(&wallB);
    // Add a data object
    data_object = new TwoWallBinForce(gflow, &wallA, &wallB);
    data_object->setMaxDistance(0.9*max_distance);
    gflow->addDataObject(data_object);
  };

  TwoWallModifier::TwoWallModifier(GFlow *gflow, const Group& group1, const Group& group2) : Modifier(gflow), wallA(WallSlideBody(gflow, group1)), 
    wallB(WallSlideBody(gflow, group2)), max_distance(0.25), acceleration(10.), dissipation(0.) {
    // Add the bodies to gflow
    gflow->addBody(&wallA);
    gflow->addBody(&wallB);
    // Add a data object
    data_object = new TwoWallBinForce(gflow, &wallA, &wallB);
    data_object->setMaxDistance(0.9*max_distance);
    gflow->addDataObject(data_object);
  };

  void TwoWallModifier::post_forces() {
    // Get wall positions.
    RealType x1 = wallA.getPosition();
    RealType x2 = wallB.getPosition();
    RealType v1 = wallA.getVelocity();
    RealType v2 = wallB.getVelocity();
    RealType dx = x2 - x1;
    gflow->minimumImage(dx, 0);
    RealType dv = v2 - v1;

    // Calculate the acceleration
    RealType a = clamp(fabs(dx) - max_distance) * (sign(dx) * acceleration - dissipation);

    // If acceleration is nonzero.
    if (a>0) {
      Vec A(sim_dimensions);
      A[0] = a;
      wallA.addAcceleration(A.data, simData);
      A.negate();
      wallB.addAcceleration(A.data, simData);
    }
  }

  void TwoWallModifier::setMaxDistance(RealType md) {
    max_distance = md;
  }

  void TwoWallModifier::setAcceleration(RealType a) {
    acceleration = a;
  }

  void TwoWallModifier::setMaxDistanceDataObject(RealType md) {
    if (data_object) data_object->setMaxDistance(md);
  }

  void TwoWallModifier::setMinDistanceDataObject(RealType md) {
    if (data_object) data_object->setMinDistance(md);
  }

  void TwoWallModifier::setBinsDataObject(int b) {
    if (data_object && b>0) data_object->setBins(b);
  }

}