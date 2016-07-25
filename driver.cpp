#include "Simulator.h"

#include <ctime>

double average(vector<vect<>> lst) {
  if (lst.empty()) return 0;
  double ave = 0;
  for (auto V : lst) ave += V.y;
  return ave/lst.size();
}

int main() {
  // Parameters
  double gap = 0.15;
  double width = 4.;
  double height = 2;
  int number = 500;
  double time = 10;
  double radius = 0.05;
  double rA = 0.05;

  // Seed random number generators
  srand48( std::time(0) );
  srand( std::time(0) );

  //----------------------------------------

  Simulator simulation;

  simulation.addStatistic(statKE);
  simulation.addStatistic(statPassiveFlow);
  simulation.addStatistic(statActiveFlow);

  simulation.setStartRecording(0);
  vect<> bias(0.1, 0);
  double act=0, pass=0;
  simulation.createControlPipe(number-25, 25, radius, 0, bias, rA, width, height);
  simulation.run(time);
  
  /// Print data
  cout << simulation.printWatchList() << endl;
  cout << "walls=" << simulation.printWalls() << ";\n";
  cout << simulation.printAnimationCommand() << endl;
  
  cout << "aveKE=" << simulation.getStatistic(0) << ";\n";
  cout << "Print[\"Average kinetic energy\"]\nListLinePlot[aveKE,PlotRange->All,PlotStyle->Black]\n";
  cout << "passflow=" << simulation.getStatistic(1) << ";\n";
  cout << "Print[\"Passive Flow\"]\nListLinePlot[passflow,PlotRange->All,PlotStyle->Black]\n";
  cout << "actflow=" << simulation.getStatistic(2) << ";\n";
  cout << "Print[\"Active Flow\"]\nListLinePlot[actflow,PlotRange->All,PlotStyle->Black]\n";
  cout << "Minimum epsilon: " << simulation.getMinEpsilon() << ";\n";
  cout << simulation.getIter() << " iterations;\n";
  cout << "Run time: " << simulation.getRunTime() << " seconds;\n";
  
  return 0;
}
