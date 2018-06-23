#include "forcemaster.hpp"
#include "force.hpp"

namespace GFlowSimulation {

  ForceMaster::ForceMaster(GFlow *gflow) : Base(gflow), ntypes(0) {};

  Force* ForceMaster::getForce(int type1, int type2) {
    if (type1>=ntypes || type2>=ntypes)
      throw ParticleTypeError("Tried types ["+toStr(type1)+"], ["+toStr(type2)+", ntypes: "+toStr(ntypes));
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
  }

}