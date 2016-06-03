#include "Simulator.h"

int main() {

  Simulator simulation;
  simulation.createHopper(100);

  simulation.run(20);

  cout << simulation.printWatchList() << endl;
  cout << "walls=" << simulation.printWalls() << ";" << endl;
  cout << "aveV=" << simulation.printAveV() << ";" << endl;
  cout << "omega=" << simulation.printNetAngularV() << ";" << endl;
  cout << "omegaSqr=" << simulation.printAveAngularVSqr() << ";" << endl;
  cout << "torque=" << simulation.printNetTorque() << ";" << endl;
  cout << "vid=" << simulation.printAnimationCommand() << endl;
  cout << "Print[\"Average Velocity\"]" << endl;
  cout << "ListLinePlot[aveV, PlotRange->All, PlotStyle->Black]\n";
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
