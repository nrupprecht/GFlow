#include "Simulator.h"

void standardHopper(Simulator &simulation) {
  simulation.discard();
  double left=0, right=1, bottom=0, top=3;
  simulation.setDimensions(left, right, bottom, top);
  double radius = 0.02;
  double gap = 0.14;
  double bottomGap = 0.05;
  double troughHeight = 0.5;
  double space = 1.0;
  double var = 0.25, mx = (1+var)*radius;
  simulation.addWall(new Wall(vect<>(0, troughHeight), vect<>(0,2*top), true));
  simulation.addWall(new Wall(vect<>(right, troughHeight), vect<>(right,2*top), true));
  simulation.addWall(new Wall(vect<>(0, troughHeight), vect<>(0.5-0.5*gap, bottomGap), true));
  simulation.addWall(new Wall(vect<>(1, troughHeight), vect<>(0.5+0.5*gap, bottomGap), true));
  simulation.addTempWall(new Wall(vect<>(0,troughHeight), vect<>(1,troughHeight), true), 3.0);
  double upper = 5; // Instead of top
  int N = 100;
  simulation.addNWParticles(N, radius, var, mx, right-mx, troughHeight+mx, upper-mx);
  simulation.setXLBound(WRAP);
  simulation.setXRBound(WRAP);
  simulation.setYTBound(NONE);
  simulation.setYBBound(RANDOM);
}

int main() {

  Simulator simulation;

  ///***** Hopper tests ***********************************************************
  srand(0);
  standardHopper(simulation);
  simulation.setSectorDims(10,10);
  simulation.run(10);
  cout << "Hopper, 100 particles, 10 seconds" << endl;
  cout << "------------------------------------------------------------\n";
  cout << "Hopper, sectors 10x10: " << simulation.getRunTime() << endl;
  cout << "Check: " << simulation.aveKE() << endl << endl;

  srand(0);
  standardHopper(simulation);
  simulation.setSectorDims(5,5);
  simulation.run(10);
  cout << "Hopper, sectors 5x5: " << simulation.getRunTime() << endl;
  cout << "Check: " << simulation.aveKE() << endl << endl;

  srand(0);
  simulation.setSectorize(false);
  simulation.run(10);
  cout << "Hopper, no sectors: " << simulation.getRunTime() << endl;
  cout << "Check: " << simulation.aveKE() << endl << endl;

  

  return 0;
}

