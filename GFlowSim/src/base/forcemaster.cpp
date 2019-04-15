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

  void ForceMaster::interact() {
    // Check if there are any interactions
    if (interactions.empty()) return;
    // Start timer
    start_timer();
    // Run all interactions
    for (auto it : interactions) it->interact();
    // Stop the timer
    stop_timer();
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

  RealType ForceMaster::getMaxCutoff(int type) const {
    return max_cutoffs[type];
  }

  RealType ForceMaster::getTotalPotentialEnergy() const {
    // Sum up all the potential energies of all the interactions.
    RealType potential = 0;
    for (auto it : interactions) potential += it->getPotential();
    // Return the total potential
    return potential;
  }

  RealType ForceMaster::getTotalVirial() const {
    // Sum up all the potential energies of all the interactions.
    RealType virial = 0;
    for (auto it : interactions) virial += it->getVirial();
    // Return the total potential
    return virial;
  }

  bool ForceMaster::typeInteracts(int type) const {
    // If type is -1, no interaction
    if (type==-1) return false;
    // Else, first make sure array is in bounds
    if (type<0 && ntypes<=type) 
      throw ParticleTypeError("From typeInteractions function.");
    // If it is, return.
    return doesInteract[type];
  }

  const vector<RealType>& ForceMaster::getMaxCutoff() const {
    return max_cutoffs;
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
    // Resize max_cutoffs array
    max_cutoffs = vector<RealType>(ntypes, 1.);
  }
  
  void ForceMaster::setInteraction(int type1, int type2, Interaction *it, bool reflexive) {
    grid[type1][type2] = it;
    if (reflexive) grid[type2][type1] = it;
    // Add to the list if it is not already there and isn't null.
    if (it!=nullptr && !contains(interactions, it)) 
      interactions.push_back(it);
    // Add to gflow's list if it is not already there and isn't null. GFlow will check for this.
    gflow->addInteraction(it);
    // Update max_cutoffs
    if (it) {
      RealType cutoff = it->getCutoff();
      if (max_cutoffs[type1]<cutoff) max_cutoffs[type1] = cutoff;
      if (max_cutoffs[type2]<cutoff) max_cutoffs[type2] = cutoff;
    }
  }

  void ForceMaster::setInteraction(Interaction *it) {
    // Get new cutoff.
    RealType cutoff = 10;
    if (it!=nullptr) { 
      // Get the cutoff
      cutoff = it->getCutoff();
      // Set interacts to true for everyone
      for (int i=0; i<ntypes; ++i) doesInteract[i] = true;
    }
    else { // Nothing interacts
      // Set interacts to false for everyone
      for (int i=0; i<ntypes; ++i) doesInteract[i] = false;
    }
    // Add to the list if it is not already there and isn't null.
    if (it!=nullptr && !contains(interactions, it)) 
      interactions.push_back(it);
    // Set the force grid.
    for (int y=0; y<ntypes; ++y) {
      for (int x=0; x<ntypes; ++x) {
        grid[y][x] = it;
      }
      // Reset the cutoff for this type.
      max_cutoffs[y] = cutoff;
    }
  }

  void ForceMaster::setCalculatePotential(bool s) {
    for (auto it : interactions) it->setDoPotential(s);
  }

  void ForceMaster::setCalculateVirial(bool s) {
    for (auto it : interactions) it->setDoVirial(s);
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
