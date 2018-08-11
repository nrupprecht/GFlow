#include "forcemaster.hpp"
// Other files
#include "force.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  ForceMaster::ForceMaster(GFlow *gflow) : Base(gflow), ntypes(0) {};

  ForceMaster::ForceMaster(GFlow *gflow, int nt) : Base(gflow) {
    setNTypes(nt);
  };

  Force* ForceMaster::getForce(int type1, int type2) {
    if (ntypes<=type1 || ntypes<=type2 || type1<0 || type2<0) return nullptr;
    return forceGrid.at(type1, type2);
  }

  void ForceMaster::clearVerletLists() {
    for (auto f : forces) f->clear();
  }

  int ForceMaster::getNTypes() const {
    return ntypes;
  }

  void ForceMaster::setNTypes(int n) {
    // Set ntypes for this object
    ntypes = n;
    // Set ntypes for the simdata
    Base::simData->ntypes = n;
    // Resize and erase array
    forceGrid.resize(n, n);
    forceGrid.setAll(nullptr);
  }

  void ForceMaster::setForce(int type1, int type2, Force *f) {
    forceGrid.at(type1, type2) = f;
    // Add to the list if it is not already there
    if (!contains(forces, f)) forces.push_back(f);
    if (!contains(gflow->forces, f)) gflow->forces.push_back(f);
  }

}
