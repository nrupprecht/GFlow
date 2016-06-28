#include "Simulator.h"

#include <ctime>

int main() {
  int number = 1;
  double time = 1;
  double radius = 0.05;

  srand48( std::time(0) );
  srand( std::time(0) );

  //----------------------------------------
  Simulator simulation;

  // simulation.setDoFluid(false);

  double delayTime = 20.;
  simulation.setMarkWatch(true);
  simulation.setDelayTime(delayTime);
  simulation.setStartTime(4.);

  // simulation.createHopper(number, radius);
  // simulation.createPipe(number);
  // simulation.createIdealGas(number, 0.03);
  // simulation.createEntropyBox(number);
  // simulation.createSquare(1);
  simulation.createFluidBox();

  simulation.setDispRate(150);
  simulation.setDispFactor(0.3);

  simulation.run(time);

  cout << "press={" << simulation.getPressurePrint() << "};\n";
  
  //cout << /* "press=" << */ simulation.printPressure() << ";\n";
  cout << "vf=" << simulation.printFV() << ";\n";

  //cout << "timeMarks=" << simulation.getTimeMarks() << ";" << endl;
  cout << simulation.printWatchList() << endl;
  cout << "walls=" << simulation.printWalls() << ";\n";
  //cout << "aveKE=" << simulation.printKE() << ";\n";
  //cout << "netOmega=" << simulation.printNetOmega() << ";\n";
  cout << simulation.printAnimationCommand() << endl;
  //cout << "Print[\"Average kinetic energy\"]\nListLinePlot[aveKE, PlotRange->All, PlotStyle->Black]\n";
  //cout << "Print[\"Net Angular Velocity\"]\nListLinePlot[netOmega, PlotRange->All, PlotStyle->Black]\n";
  //cout << simulation.getMinEpsilon() << ";\n";
  //cout << simulation.getIter() << ";\n";
  //cout << simulation.getRunTime() << ";";
  
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

