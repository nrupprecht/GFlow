#include "bond.hpp"

namespace GFlowSimulation {

  Bond::Bond(GFlow *gflow) : Modifier(gflow) {};

  void Bond::addBond(int id1, int id2) {
    left.push_back(id1);
    right.push_back(id2);
  }

}