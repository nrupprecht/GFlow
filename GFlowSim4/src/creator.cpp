#include "creator.hpp"

namespace GFlowSimulation {

  Creator::Creator(int ac, char **av) : argc(ac), argv(av), ourParser(true) {
    parserPtr = new ArgParse(ac, av);
  };

  Creator::Creator(ArgParse *p) : argc(p->getArgc()), argv(p->getArgv()), parserPtr(p), ourParser(false) {};

  Creator::~Creator() {
    if (ourParser && parserPtr) delete parserPtr;
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

  void Creator::hs_relax(GFlow* gflow, RealType time) {
    // Check for valid object
    if (gflow==nullptr) return;

    // Use overdamped integrator for relaxation
    Integrator *integrator = gflow->integrator; // Save integrator
    gflow->integrator = new OverdampedIntegrator(gflow);

    // Use hard sphere forces
    int ntypes = gflow->getNTypes();
    HardSphere hsForce(gflow);
    ForceMaster *master = gflow->forceMaster; // Save old master
    vector<Interaction*> interactions = gflow->interactions; // Save old forces
    gflow->interactions.clear(); // Clear old forces
    // New force master - for hard sphere forces only.
    ForceMaster forceMaster(gflow, ntypes);
    gflow->forceMaster = &forceMaster; // Give it to gflow
    // All particles interact as hard spheres
    for (int n1=0; n1<ntypes; ++n1) 
      for (int n2=0; n2<ntypes; ++n2)
        forceMaster.setForce(n1, n2, &hsForce);

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

    // Reset forces
    gflow->forceMaster = master;
    gflow->interactions = interactions;
    gflow->simData->clearF();
  }

  void Creator::relax(class GFlow *gflow, RealType time) {
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
