#include "bond.hpp"

namespace GFlowSimulation {

  Bond::Bond(GFlow *gflow) : Modifier(gflow) {};

  void Bond::pre_integrate() {
    updateLocalIDs();
  }

  bool Bond::checkBondVectors() const {
    return left.size()==right.size() && left.size()==gleft.size() && gleft.size()==gright.size();
  }

  void Bond::updateLocalIDs() {
    // Make sure sizes are the same
    if (!checkBondVectors()) throw UnequalBondVectors();
    int nbonds = left.size();
    // Update local ids
    SimData *sd = Base::simData;
    for (int i=0; i<nbonds; ++i) {
      int gid1 = gleft[i], gid2 = gright[i];
      int id1 = sd->getLocalID(gid1), id2 = sd->getLocalID(gid2);
      left[i] = id1; 
      right[i] = id2;
    }
  }

}