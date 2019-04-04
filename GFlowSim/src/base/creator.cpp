#include "creator.hpp"
// Other files
#include "../interactions/interaction-choice.hpp"

namespace GFlowSimulation {

  Creator::Creator(int ac, char **av) : argc(ac), argv(av), ourParser(true), simBounds(Bounds(2)), sim_dimensions(2) {
    parserPtr = new ArgParse(ac, av);
  };

  Creator::Creator(ArgParse *p) : argc(p->getArgc()), argv(p->getArgv()), parserPtr(p), ourParser(false), simBounds(Bounds(2)), sim_dimensions(2) {};

  Creator::~Creator() {
    if (ourParser && parserPtr) 
      if (parserPtr) delete parserPtr;
  }

  void Creator::setSeed(uint s) {
    seed = s;
  }

  unsigned Creator::getSeed() {
    return seed;
  }

  void Creator::seedGenerator(uint s) {
    seed = s;
  };

  void Creator::setBCFlag(BCFlag b) { 
    bcFlag = b; 
  }

  void Creator::setDimensions(int d) {
    sim_dimensions = d;
  }
  
  void Creator::hs_relax(GFlow* gflow, RealType time, bool relax_integrator) {
    // Check for valid object
    if (gflow==nullptr) return;
    
    // Get sim dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Save gflow's data
    Integrator *integrator = gflow->integrator; // Save integrator
    ForceMaster *master    = gflow->forceMaster; // Save old master
    vector<Interaction*> interactions = gflow->interactions; // Save old forces

    // Use hard sphere forces
    int ntypes = gflow->getNTypes();
    Interaction *hsForce = InteractionChoice::choose(gflow, InteractionChoice::HardSphereToken, sim_dimensions);
    
    // New force master - has only hard sphere forces.
    ForceMaster forceMaster(gflow, ntypes);
    Integrator *rx_integrator;
    if (relax_integrator) rx_integrator = new RelaxIntegrator(gflow);
    else rx_integrator = new OverdampedIntegrator(gflow);

    // Give gflow new data
    gflow->interactions.clear(); // Clear old forces
    gflow->integrator = rx_integrator;
    gflow->forceMaster = &forceMaster; // Give it to gflow

    // All particles interact as hard spheres
    for (int n1=0; n1<ntypes; ++n1) 
      for (int n2=0; n2<ntypes; ++n2)
        forceMaster.setInteraction(n1, n2, hsForce);

    // Make sure all forces are zero
    gflow->simData->clearF();
    // Relax simulation
    gflow->requestTime(time);
    gflow->run(); // GFlow run calls initialize, so all the base objects' pointers will be correct.
    // Reset times
    gflow->dataMaster->resetTimer(); // So the setup does not count towards the run time.
    gflow->resetAllTimes();

    // Reset GFlow
    gflow->integrator   = integrator;
    gflow->forceMaster  = master;
    gflow->interactions = interactions;
    gflow->initialize(); // Correct all the base objects' pointers.

    // Clear forces
    gflow->simData->clearF();

    // Clean up
    delete rx_integrator;
    delete hsForce;
  }

  void Creator::relax(class GFlow *gflow, RealType time) {
    // Check for valid object
    if (gflow==nullptr) return;

    // Use overdamped integrator for relaxation
    Integrator *integrator = gflow->integrator; // Save integrator
    RelaxIntegrator rx_integrator(gflow);
    gflow->integrator = &rx_integrator;

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Relax simulation
    gflow->requestTime(time); // Should be long enough
    gflow->run();
    gflow->dataMaster->resetTimer(); // So the setup does not count towards the run time
    gflow->resetAllTimes();
    
    // Reset integrator
    gflow->integrator = integrator;
  }

  void Creator::fix_particle_velocities(SimData *simData) {
    for (auto &fix : particle_fixers) {
      // Get local id of the particle
      int id = simData->getLocalID(fix.global_id);
      if (0<=id && id<simData->size()) {
        // Set the initial velocity of the particle.
        copyVec(fix.velocity, simData->V(id), sim_dimensions);
      }
    }
  }

  void Creator::clear_particle_fixers() {
    particle_fixers.clear();
  }

}
