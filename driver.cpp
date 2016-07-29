#include "Simulator.h"

#include <ctime>

int main() {
  // Parameters
  double width = 4.;
  double height = 2.;
  double Vol = width*height;
  int number = 900-45; // 900 -> max
  double time = 10;
  double radius = 0.05;
  double rA = radius;
  double pA = 0.1;
  int NA = number*pA, NP = number-NA;

  // Seed random number generators
  srand48( std::time(0) );
  srand( std::time(0) );

  //----------------------------------------

  Simulator simulation;

  simulation.addStatistic(statKE);
  simulation.addStatistic(statPassiveFlow);
  simulation.addStatistic(statActiveFlow);
  simulation.addStatistic(statFlowRatio);

  simulation.setStartRecording(3);
  vect<> bias(0.1, 0);
  double act=0, pass=0;
  simulation.createControlPipe(NP, NA, radius, 0, bias, rA, width, height);
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
  cout << "xi=" << simulation.getStatistic(3) << ";\n";
  cout << "Print[\"Xi\"]\nListLinePlot[xi,PlotRange->All,PlotStyle->Black]\n";

  cout << "Minimum epsilon: " << simulation.getMinEpsilon() << ";\n";
  cout << simulation.getIter() << " iterations;\n";
  cout << "Run time: " << simulation.getRunTime() << " seconds;\n";
  
  cout << "Packing ratio: " << number*PI*sqr(radius)/(width*height) << endl;

  return 0;
}
