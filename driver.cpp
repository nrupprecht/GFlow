#include "Simulator.h"

#include <ctime>

int main() {
  // Parameters
  double width = 4.; // 4
  double height = 4.; // 2
  double radius = 0.05;
  double velocity = 0.5;
  double phi = 0.5;
  double time = 300;
  double start = 0;
  double pA = 0.;

  // Display parameters
  bool animate = false;
  bool dispKE = false;
  bool dispFlow = false;
  bool dispProfile = false;
  bool dispAveProfile = true;

  // Dependent variables
  double Vol = width*height;
  double rA = radius;
  int number = Vol/(PI*sqr(radius))*phi;
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

  simulation.setStartRecording(start);
  vect<> bias(0.1, 0);
  double act=0, pass=0;
  simulation.createControlPipe(NP, NA, radius, velocity, bias, rA, width, height);
  simulation.run(time);

  /// Print condition summary
  cout << "Dimensions: " << width << " x " << height << "\n";
  cout << "Radius: " << radius << "\n";
  cout << "Fluid Velocity: " << velocity << "\n";
  cout << "Phi: " << phi << ", Number: " << number << ", (Actual Phi: " << number*PI*sqr(radius)/(width*height) << ")\n";
  cout << "Percent Active: " << pA*100 << "%, (Actual %: " << 100.*NA/(double)(NA+NP) << ")\n";
  cout << "N Active: " << NA << ", N Passive: " << NP << "\n";
  cout << "Sim Time: " << time << ", Run time: " << simulation.getRunTime() << ", Ratio: " << time/simulation.getRunTime() << endl;
  cout << "Start Time: " << start << "\n";
  cout << "Iters: " << simulation.getIter() << "\n\n";
  

  /// Print data
  if (animate) {
    cout << simulation.printWatchList() << endl;
    cout << "walls=" << simulation.printWalls() << ";\n";
    cout << simulation.printAnimationCommand() << endl;
  }
  
  if (dispKE) {
    cout << "aveKE=" << simulation.getStatistic(0) << ";\n";
    cout << "Print[\"Average kinetic energy\"]\nListLinePlot[aveKE,PlotRange->All,PlotStyle->Black]\n";
  }
  if (dispFlow) {
    cout << "passflow=" << simulation.getStatistic(1) << ";\n";
    cout << "Print[\"Passive Flow\"]\nListLinePlot[passflow,PlotRange->All,PlotStyle->Black]\n";
    cout << "actflow=" << simulation.getStatistic(2) << ";\n";
    cout << "Print[\"Active Flow\"]\nListLinePlot[actflow,PlotRange->All,PlotStyle->Black]\n";
    cout << "xi=" << simulation.getStatistic(3) << ";\n";
    cout << "Print[\"Xi\"]\nListLinePlot[xi,PlotRange->All,PlotStyle->Black]\n";
  }
  if (dispProfile) {
    cout << "prof=" << simulation.getProfile() << ";\n";
    cout << "Print[\"Profile\"]\nMatrixPlot[prof,ImageSize->Large]\n"; 
  }
  if (dispAveProfile) {
    cout << "aveProf=" << simulation.getAveProfile() << ";\n";
    cout << "Print[\"Average Profile\"]\nListLinePlot[aveProf,PlotStyle->Black]\n";
  }

  return 0;
}
