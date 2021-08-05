#include <bonded/twogroupharmonic.hpp>

using namespace GFlowSimulation;

TwoGroupHarmonic::TwoGroupHarmonic(GFlow *gflow)
    : Bonded(gflow), groupA(gflow), groupB(gflow), springConstant(50.),
      min_distance(0.), max_distance(0.15) {};

TwoGroupHarmonic::TwoGroupHarmonic(GFlow *gflow, Group &gA, Group &gB)
    : Bonded(gflow), groupA(gA), groupB(gB),
      springConstant(50.), min_distance(0.),
      max_distance(0.15) {};

void TwoGroupHarmonic::pre_integrate() {
  groupA.update_local_ids();
  groupB.update_local_ids();

  massA = groupA.findTotalMass();
  massB = groupB.findTotalMass();
}

void TwoGroupHarmonic::interact() const {
  if (groupA.empty() || groupB.empty() || massA <= 0 || massB <= 0) {
    return;
  }

  if (simData->getNeedsRemake()) {
    groupA.update_local_ids();
    groupB.update_local_ids();
  }

  // Compute the distance between the centers of mass of the groups. Use this to calculate the force between them.
  Vec X1(sim_dimensions), X2(sim_dimensions), dX(sim_dimensions), F(sim_dimensions);
  groupA.findCenterOfMass(X1.data);
  groupB.findCenterOfMass(X2.data);
  gflow->getDisplacement(X2.data, X1.data, dX.data);
  RealType d = magnitude(dX);
  dX.normalize();
  // Compute accelerations
  if (d < min_distance) {
    F = springConstant * (d - min_distance) * dX;
    groupA.addAcceleration((1. / massA * F).data);
    F.negate();
    groupB.addAcceleration((1. / massB * F).data);
  }
  if (max_distance < d) {
    F = springConstant * (d - max_distance) * dX;
    groupA.addAcceleration((1. / massA * F).data);
    F.negate();
    groupB.addAcceleration((1. / massB * F).data);
  }
}

void TwoGroupHarmonic::setMaxDistance(RealType md) {
  max_distance = md;
}

void TwoGroupHarmonic::setMinDistance(RealType md) {
  min_distance = md;
}

void TwoGroupHarmonic::setSpringConstant(RealType a) {
  springConstant = a;
}

void TwoGroupHarmonic::setGroupA(Group &gA) {
  groupA = gA;
  massA = groupA.findTotalMass();
}

void TwoGroupHarmonic::setGroupB(Group &gB) {
  groupB = gB;
  massB = groupA.findTotalMass();
}

int TwoGroupHarmonic::size() const {
  return groupA.size() + groupB.size();
}
