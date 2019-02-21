#include "angle.hpp"

namespace GFlowSimulation {

  Angle::Angle(GFlow *gflow) : Modifier(gflow) {};

  void Angle::addAngle(int id1, int id2, int id3) {
    left.push_back(id1);
    center.push_back(id2);
    right.push_back(id3);
  }

}