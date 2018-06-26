#include "creator.hpp"

namespace GFlowSimulation {

  Creator::Creator(int ac, char **av) : argc(ac), argv(av), ourParser(true) {
    parserPtr = new ArgParse(ac, av);
  };

  Creator::Creator(ArgParse *p) : argc(p->getArgc()), argv(p->getArgv()), parserPtr(p), ourParser(false) {};

  Creator::~Creator() {
    if (ourParser && parserPtr) delete parserPtr;
  }
}