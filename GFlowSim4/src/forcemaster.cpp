#include "forcemaster.hpp"
#include "force.hpp"

namespace GFlowSimulation {

  ForceMaster::ForceMaster(GFlow *gflow) : Base(gflow), ntypes(0) {};

  Force* ForceMaster::getForce(int type1, int type2) {
    if (ntypes<=type1 || ntypes<=type2 || type1<0 || type2<0) return nullptr;
    return forceGrid.at(type1, type2);
  }

  void ForceMaster::clearVerletLists() {
    for (auto f : forces) f->clearVerletList();
  }

  void ForceMaster::setNTypes(int n) {
    ntypes = n;
    forceGrid.resize(n, n);
  }

  void ForceMaster::setForce(int type1, int type2, Force *f) {
    forceGrid.at(type1, type2) = f;
    // Add to the list if it is not already there
    if (!contains(forces, f)) forces.push_back(f);
    if (!contains(gflow->forces, f)) gflow->forces.push_back(f);
  }

}
