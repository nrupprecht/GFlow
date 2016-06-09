#include "Simulator.h"

int main() {
  int number = 100;
  double time = 10;
  double radius = 0.03;

  //----------------------------------------
  Simulator simulation;
  simulation.createHopper(number, radius);
  // simulation.createPipe(100);
  // simulation.createIdealGas(100, 0.0075);
  // simulation.createEntropyBox(100);
  // simulation.createSquare(1);

  simulation.setSectorize(true);
  simulation.setAdjustEpsilon(true);
  //simulation.setMinEpsilon(1e-4);

  simulation.setDispRate(50);
  //simulation.setDispFactor(0.3);
  simulation.run(time);
  
  cout << simulation.printWatchList() << endl;
  cout << "walls=" << simulation.printWalls() << ";\n";
  cout << "aveKE=" << simulation.printKE() << ";\n";
  cout << "omegaSqr=" << simulation.printAveOmegaSqr() << ";\n";
  cout << simulation.printAnimationCommand() << endl;
  cout << "Print[\"Average kinetic energy\"]\n";
  cout << "ListLinePlot[aveKE, PlotRange->All, PlotStyle->Black]\n";
  cout << "Print[\"Average Angular Velocity Squared\"]\n";
  cout << "ListLinePlot[omegaSqr, PlotRange->All, PlotStyle->Black]\n";
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

