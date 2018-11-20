#include "creator.hpp"

namespace GFlowSimulation {

  Creator::Creator(int ac, char **av) : argc(ac), argv(av), ourParser(true) {
    parserPtr = new ArgParse(ac, av);
  };

  Creator::Creator(ArgParse *p) : argc(p->getArgc()), argv(p->getArgv()), parserPtr(p), ourParser(false) {};

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

  void Creator::hs_relax(GFlow* gflow, RealType time) const {
    // Check for valid object
    if (gflow==nullptr) return;

    // Save gflow's data
    Integrator *integrator = gflow->integrator; // Save integrator
    ForceMaster *master = gflow->forceMaster; // Save old master
    vector<Interaction*> interactions = gflow->interactions; // Save old forces

    // Use hard sphere forces
    int ntypes = gflow->getNTypes();
    HardSphere hsForce(gflow);
    // New force master - has only hard sphere forces.
    ForceMaster forceMaster(gflow, ntypes);
    
    // Give gflow new data
    gflow->interactions.clear(); // Clear old forces
    gflow->integrator = new OverdampedIntegrator(gflow);
    gflow->forceMaster = &forceMaster; // Give it to gflow

    // All particles interact as hard spheres
    for (int n1=0; n1<ntypes; ++n1) 
      for (int n2=0; n2<ntypes; ++n2)
        forceMaster.setInteraction(n1, n2, &hsForce);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Relax simulation
    gflow->requestTime(time);
    gflow->run();
    // Reset times
    gflow->dataMaster->resetTimer(); // So the setup does not count towards the run time
    gflow->resetAllTimes();

    // Reset overdamped integrator
    delete gflow->integrator;

    // Reset GFlow
    gflow->integrator = integrator;
    gflow->forceMaster = master;
    gflow->interactions = interactions;

    // Clear forces
    gflow->simData->clearF();
  }

  void Creator::relax(class GFlow *gflow, RealType time) const {
    // Check for valid object
    if (gflow==nullptr) return;

    // Use overdamped integrator for relaxation
    Integrator *integrator = gflow->integrator; // Save integrator
    gflow->integrator = new OverdampedIntegrator(gflow);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Relax simulation
    gflow->requestTime(time); // Should be long enough
    gflow->run();
    gflow->dataMaster->resetTimer(); // So the setup does not count towards the run time
    gflow->resetAllTimes();
    
    // Reset integrator
    delete gflow->integrator;
    gflow->integrator = integrator;

  }

}
