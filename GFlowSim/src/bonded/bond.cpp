#include "bond.hpp"

namespace GFlowSimulation {

  Bond::Bond(GFlow *gflow) : Bonded(gflow) {};

  void Bond::pre_integrate() {
    updateLocalIDs();
  }

  int Bond::size() const {
    return left.size();
  }

  bool Bond::checkBondVectors() const {
    return left.size()==right.size() && left.size()==gleft.size() && gleft.size()==gright.size();
  }

  void Bond::updateLocalIDs() const {
    // Make sure sizes are the same
    if (!checkBondVectors()) throw UnequalBondVectors();
    int nbonds = left.size();
    // Update local ids
    for (int i=0; i<nbonds; ++i) {
      int gid1 = gleft[i], gid2 = gright[i];
      int id1 = simData->getLocalID(gid1), id2 = simData->getLocalID(gid2);
      left[i] = id1; 
      right[i] = id2;
    }
  }

}