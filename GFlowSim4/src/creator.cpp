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

  void Creator::relax(GFlow* gflow, RealType time) {
    // Check for valid object
    if (gflow==nullptr) return;

    // Use overdamped integrator for relaxation
    gflow->integrator = new OverdampedIntegrator(gflow);

    // Make sure all forces are zero
    gflow->simData->clearF();

    // Relax simulation
    gflow->requestTime(time); // Should be long enough
    gflow->run();
    gflow->dataMaster->resetTimer(); // So the setup does not count towards the run time

    // Set new integrator
    delete [] gflow->integrator;

    // Make sure all forces are zero again
    gflow->simData->clearF();
  }

}