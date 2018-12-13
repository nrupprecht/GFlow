#include "forcemaster.hpp"
// Other files
#include "interaction.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  ForceMaster::ForceMaster(GFlow *gflow) : Base(gflow), ntypes(0) {};

  ForceMaster::ForceMaster(GFlow *gflow, int nt) : Base(gflow) {
    setNTypes(nt);
  };

  Interaction* ForceMaster::getInteraction(int type1, int type2) {
    if (ntypes<=type1 || ntypes<=type2 || type1<0 || type2<0) return nullptr;
    //return grid.at(type1, type2);
    return grid[type1][type2];
  }

  void ForceMaster::clear() {
    for (auto it : interactions) it->clear();
  }

  void ForceMaster::close() {
    for (auto it : interactions) it->close();
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
    grid = vector<vector<Interaction*> >(n, vector<Interaction*>(n, nullptr));
  }

  void ForceMaster::setInteraction(int type1, int type2, Interaction *it) {
    //grid.at(type1, type2) = it;
    grid[type1][type2] = it;
    // Add to the list if it is not already there
    if (!contains(interactions, it)) 
      interactions.push_back(it);
    // Add to gflow's list if it is not already there
    if (!contains(gflow->interactions, it)) 
      gflow->interactions.push_back(it);
  }

}
