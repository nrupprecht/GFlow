#include "GFlowBase.h"

int main () {

  cout << sizeof(Particle) << endl;
  cout << sizeof(Wall) << endl;

  GFlowBase simulator;

  simulator.addParticle(0,0,1);

  simulator.run(1.);

  return 0;
}
