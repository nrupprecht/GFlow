#include "Simulator.h"

int main() {
  int number = 1000;
  double time = 20;
  double radius = 0.03;

  //----------------------------------------
  Simulator simulation;
  simulation.createHopper(number, radius);
  // simulation.createPipe(number);
  // simulation.createIdealGas(number, 0.03);
  // simulation.createEntropyBox(number);
  // simulation.createSquare(1);

  simulation.setDispRate(50);
  //simulation.setDispFactor(0.3);
  simulation.run(time);
  
  cout << simulation.printWatchList() << endl;
  cout << "walls=" << simulation.printWalls() << ";\n";
  cout << "aveKE=" << simulation.printKE() << ";\n";
  //cout << "netOmega=" << simulation.printNetOmega() << ";\n";
  cout << simulation.printAnimationCommand() << endl;
  cout << "Print[\"Average kinetic energy\"]\nListLinePlot[aveKE, PlotRange->All, PlotStyle->Black]\n";
  //cout << "Print[\"Net Angular Velocity\"]\nListLinePlot[netOmega, PlotRange->All, PlotStyle->Black]\n";
  cout << simulation.getMinEpsilon() << ";\n";
  cout << simulation.getIter() << ";\n";
  cout << simulation.getRunTime() << ";";
  
  return 0;
}

// COMPARISON: 10 SECOND SIMULATION, IDEAL GAS

// S , 100 - 39s
// NS, 100 - 132s

// S , 200 - 111s
// NS, 200 - 404s

// S , 250 - 140s
// NS, 250 - 572s

// S , 500 - 336s
// NS, 500 - 1509s

