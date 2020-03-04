#include "forcemaster.hpp"
// Other files
#include "interaction.hpp"
#include "integrator.hpp"
#include "simdata.hpp"
#include "topology.hpp"

namespace GFlowSimulation {

  ForceMaster::ForceMaster(GFlow *gflow) : Base(gflow) {
    // One type by default
    setNTypes(1);
  }

  ForceMaster::ForceMaster(GFlow *gflow, int nt) : Base(gflow) {
    setNTypes(nt);
  };

  void ForceMaster::initialize() {
    // Base initialization.
    Base::initialize();
    // Initialize does interact.
    initialize_does_interact();
  }

  void ForceMaster::pre_integrate() {
    // Clear timer.
    clear_timer();
    // Initialize does interect. \todo Probably unneccesary, already done in initialize.
    initialize_does_interact();
    // Set Max DT if it is smaller than the preexisting max_dt.
    if (integrator) {
      compute_timescale();
      RealType suggested_dt = time_scale*time_scale_factor;
      if (0<suggested_dt && suggested_dt<integrator->getMaxDT()) integrator->setMaxDT(suggested_dt);
    }
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

  void ForceMaster::interact_ghosts() {
    // Check if there are any interactions
    if (interactions.empty() || topology->getNumProc()==1) return;
    // Start timer
    start_timer();
    // Run all interactions
    for (auto it : interactions) it->interact_ghosts();
    // Stop the timer
    stop_timer();
  }

  shared_ptr<Interaction> ForceMaster::getInteraction(int type1, int type2) {
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
    // Return the total potential.
    return virial;
  }

  bool ForceMaster::interactingTypes() const {
    return any_interactions;
  }

  bool ForceMaster::typeInteracts(int type) const {
    // If type is -1, no interaction.
    if (type==-1) return false;
    // Else, first make sure array is in bounds.
    if (type<0 || ntypes<=type) 
      throw ParticleTypeError("From ForceMaster::typeInteracts(int) function. Type=" +toStr(type) + ". ntypes=" + toStr(ntypes));
    // If it is, return.
    return doesInteract[type];
  }

  bool ForceMaster::typesInteract(int type1, int type2) const {
    // If type is -1, no interaction
    if (type1==-1 || type2==-1) return false;
    // Else, first make sure array is in bounds
    if (type1<0 || ntypes<=type1 || type2<0 || ntypes<type2) 
      throw ParticleTypeError("From ForceMaster::typeInteracts(int, int) function. Type1=" + toStr(type1) + ", Type2=" + toStr(type2) + ". ntypes=" + toStr(ntypes));
    // If it is, return.
    return grid[type1][type2]!=nullptr;
  }

  const vector<RealType>& ForceMaster::getMaxCutoff() const {
    return max_cutoffs;
  }

  void ForceMaster::setNTypes(int n) {
    // Set ntypes for this object
    ntypes = n;
    if (simData) simData->_ntypes = n;
    // Resize and erase array. Use decltype since the type of grid is very long.
    grid = decltype(grid)(n, vector<shared_ptr<Interaction> >(n, nullptr));
    // Resize and reset doesInteract
    doesInteract = vector<bool>(ntypes, false);
    for (int i=0; i<ntypes; ++i) doesInteract[i] = false; // Since all interactions in the grid are null
    // Resize max_cutoffs array
    max_cutoffs = vector<RealType>(ntypes, 1.);
  }
  
  void ForceMaster::setInteraction(int type1, int type2, shared_ptr<Interaction> it, bool reflexive) {
    // Set interaction between type1 and type2 to be *it.
    grid[type1][type2] = it;
    // If reflexive, then the interaction doesn't care which type of particle is "first" in the pairing,
    // and type2 interacts with type1 via the same interaction.
    if (reflexive) grid[type2][type1] = it;
    // Add to the list if it is not already there and isn't null.
    if (it && !contains(interactions, it)) interactions.push_back(it);
    // Add to gflow's list if it is not already there and isn't null. GFlow will check for this.
    gflow->addInteraction(it);
    // Update max_cutoffs
    if (it) {
      RealType cutoff = it->getCutoff();
      if (max_cutoffs[type1]<cutoff) max_cutoffs[type1] = cutoff;
      if (max_cutoffs[type2]<cutoff) max_cutoffs[type2] = cutoff;
    }
  }

  void ForceMaster::setInteraction(shared_ptr<Interaction> it) {
    for (int i=0; i<ntypes; ++i)
      for (int j=i; j<ntypes; ++j)
        setInteraction(i, j, it, true);
  }

  void ForceMaster::setCalculatePotential(bool s) {
    for (auto it : interactions) it->setDoPotential(s);
  }

  void ForceMaster::setCalculateVirial(bool s) {
    for (auto it : interactions) it->setDoVirial(s);
  }

  void ForceMaster::initialize_does_interact() {
    // Set up doesInteract array - default value is false.
    for (int i=0; i<ntypes; ++i) 
      doesInteract[i] = false;
    // Check which values of doesInteract should be true.
    for (int i=0; i<ntypes; ++i) 
      for (int j=0; j<ntypes; ++j) {
        if (grid[i][j]!=nullptr) {
          // Types i and j both interact with at least one other type.
          doesInteract[i]  = true;
          doesInteract[j]  = true;
          // There is at least one interaction.
          any_interactions = true;
        }
      }
  }

  inline void ForceMaster::compute_timescale() {
    // Need valid objects.
    if (simData==nullptr) return;
    // Find minimum mass particle.
    RealType im = 0;
    for (int i=0; i<simData->size_owned(); ++i)
      if (simData->Im(i)>im) im = simData->Im(i);
    if (im==0) return;
    // Set min mass.
    RealType mass = 1./im;

    // Reset time scale
    time_scale = -1.;
    RealType suggestion = 0;
    // Check all interactions for suggestions.
    for (auto it : interactions) {
      if (it) {
        suggestion = it->suggest_timescale(mass);
        if (0<suggestion && (suggestion<time_scale || time_scale==-1)) time_scale = suggestion;
      }
    }
  }

}
