#include "bond.hpp"

namespace GFlowSimulation {

  Bond::Bond(GFlow *gflow) : Modifier(gflow) {};

  bool Bond::checkBondVectors() {
    return left.size()==right.size() && left.size()==gleft.size() && gleft.size()==gright.size();
  }

}