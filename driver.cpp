#include "Simulator.h"

int main() {

  Simulator simulation;
  simulation.createHopper();

  simulation.run(10);

  cout << simulation.printWatchList() << endl;
  cout << "walls=" << simulation.printWalls() << ";" << endl;
  cout << "aveV=" << simulation.printAveV() <<";" << endl;
  cout << "maxV=" << simulation.printMaxV() << ";" << endl;
  cout << "vid=" << simulation.printAnimationCommand() << endl;
  cout << "Print[\"Average Velocity\"]" << endl;
  cout << "ListLinePlot[aveV, PlotRange->All, PlotStyle->Black]\n";
  cout << "Print[\"Maximum Velocity\"]" << endl;
  cout << "ListLinePlot[maxV, PlotRange->All, PlotStyle->Black]\n";
  cout << simulation.getMinEpsilon() << ";\n";
  cout << simulation.getIter() << ";";
  
  return 0;
}
