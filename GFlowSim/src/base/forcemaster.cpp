#include "forcemaster.hpp"
// Other files
#include "interaction.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  ForceMaster::ForceMaster(GFlow *gflow) : Base(gflow) {
    // One type by default
    setNTypes(1);
  }

  ForceMaster::ForceMaster(GFlow *gflow, int nt) : Base(gflow) {
    setNTypes(nt);
  };

  void ForceMaster::initialize() {
    // Base initialization
    Base::initialize();
    // Initialize
    initialize_does_interact();
  }

  void ForceMaster::pre_integrate() {
    initialize_does_interact();
  }

  Interaction* ForceMaster::getInteraction(int type1, int type2) {
    if (ntypes<=type1 || ntypes<=type2 || type1<0 || type2<0) return nullptr;
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
    if (simData) simData->_ntypes = n;
    // Resize and erase array
    grid = vector<vector<Interaction*> >(n, vector<Interaction*>(n, nullptr));
    // Resize and reset doesInteract
    if (doesInteract) delete [] doesInteract;
    doesInteract = new bool[ntypes];
    for (int i=0; i<ntypes; ++i) doesInteract[i] = false; // Since all interactions in the grid are null
  }
  
  void ForceMaster::setInteraction(int type1, int type2, Interaction *it) {
    //grid.at(type1, type2) = it;
    grid[type1][type2] = it;
    // Add to the list if it is not already there and isn't null.
    if (it!=nullptr && !contains(interactions, it)) 
      interactions.push_back(it);
    // Add to gflow's list if it is not already there and isn't null. GFlow will check for this.
    gflow->addInteraction(it);
  }

  bool ForceMaster::typeInteracts(int type) {
    // If type is -1, no interaction
    if (type==-1) return false;
    // Else, first make sure array is in bounds
    if (type<0 && ntypes<=type) 
      throw ParticleTypeError("From typeInteractions function.");
    // If it is, return.
    return doesInteract[type];
  }

  void ForceMaster::initialize_does_interact() {
    // Set up doesInteract array - default value is false
    for (int i=0; i<ntypes; ++i) 
      doesInteract[i] = false;
    // Check which values of doesInteract should be true
    for (int i=0; i<ntypes; ++i) 
      for (int j=0; j<ntypes; ++j) {
        if (grid[i][j]!=nullptr) {
          doesInteract[i] = true;
          doesInteract[j] = true;
        }
      }
  }

}
