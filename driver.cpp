#include "Simulator.h"

int main() {

  Simulator simulation;
  //simulation.createHopper(250);
  //simulation.createPipe(100);
  simulation.createIdealGas(300);

  simulation.setDispRate(60);
  //simulation.setDispFactor(0.3);
  simulation.run(30);

  cout << simulation.printWatchList() << endl;
  cout << "walls=" << simulation.printWalls() << ";" << endl;
  cout << "aveVSqr=" << simulation.printAveVSqr() << ";" << endl;
  cout << "omega=" << simulation.printNetOmega() << ";" << endl;
  cout << "omegaSqr=" << simulation.printAveOmegaSqr() << ";" << endl;
  cout << "torque=" << simulation.printNetTorque() << ";" << endl;
  cout << simulation.printAnimationCommand() << endl;
  cout << "Print[\"Average Velocity Squared\"]" << endl;
  cout << "ListLinePlot[aveVSqr, PlotRange->All, PlotStyle->Black]\n";
  cout << "Print[\"Angular Velocity\"]" << endl;
  cout << "ListLinePlot[omega, PlotRange->All, PlotStyle->Black]\n";
  cout << "Print[\"Average Angular Velocity Squared\"]" << endl;
  cout << "ListLinePlot[omegaSqr, PlotRange->All, PlotStyle->Black]\n";
  cout << "Print[\"Torque\"]" << endl;
  cout << "ListLinePlot[torque, PlotRange->All, PlotStyle->Black]\n";
  cout << simulation.getMinEpsilon() << ";\n";
  cout << simulation.getIter() << ";";
  
  return 0;
}
