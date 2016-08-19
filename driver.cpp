#include "Simulator.h"

int main(int argc, char** argv) {
  auto start_t = clock();
  // Parameters
  double width = 2.; // 4
  double height = 1.; // 2
  double radius = 0.05;
  double velocity = 0.5;
  double phi = 0.5;
  double time = 600.;
  double start = 30;
  double pA = 0.;
  double activeF = 0.1; // default is 5
  int samplePoints = -1; 

  // Display parameters
  bool animate = false;
  bool dispKE = false;
  bool aveKE = false;
  bool dispFlow = false;
  bool dispProfile = false;
  bool dispAveProfile = true;

  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser(argc, argv);
  pair<string,string> opt;
  stringstream stream;
  opt = parser.find("width");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> width;
  }
  opt = parser.find("height");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> height;
  }
  opt = parser.find("radius");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> radius;
  }
  opt = parser.find("velocity");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> velocity;
  }
  opt = parser.find("phi");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> phi;
  }
  opt = parser.find("time");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> time;
  }
  opt = parser.find("start");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> start;
  }
  opt = parser.find("active");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> pA;
  }
  opt = parser.find("force");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> activeF;
  }
  opt = parser.find("animate");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> animate;
  }
  opt = parser.find("dispKE");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> dispKE;
  }
  opt = parser.find("aveKE");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> aveKE;
  }
  opt = parser.find("flow");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> dispFlow;
  }
  opt = parser.find("profileMap");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> dispProfile;
  }
  opt = parser.find("profile");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> dispAveProfile;
  }
  opt = parser.find("points");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> samplePoints;
  }
  //----------------------------------------

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
  simulation.addStatistic(statPassiveKE);
  simulation.addStatistic(statPassiveFlow);
  simulation.addStatistic(statActiveFlow);
  simulation.addStatistic(statFlowRatio);
  simulation.setStartRecording(start);
  simulation.createControlPipe(NP, NA, radius, velocity, activeF, rA, width, height);
  if (samplePoints>0) simulation.setSamplePoints(samplePoints);
  simulation.run(time);
  auto end_t = clock();
  
  /// Print condition summary
  cout << "Dimensions: " << width << " x " << height << "\n";
  cout << "Radius: " << radius << "\n";
  cout << "Fluid Velocity: " << velocity << "\n";
  cout << "Phi: " << phi << ", Number: " << number << ", (Actual Phi: " << number*PI*sqr(radius)/(width*height) << ")\n";
  cout << "Percent Active: " << pA*100 << "%, (Actual %: " << 100.*NA/(double)(NA+NP) << ")\n";
  if (NA>0) cout << "Active Force: " << activeF << endl;
  cout << "N Active: " << NA << ", N Passive: " << NP << "\n";
  cout << "Sim Time: " << time << ", Run time: " << simulation.getRunTime() << ", Ratio: " << time/simulation.getRunTime() << endl;
  cout << "Start Time: " << start << "\n";
  cout << "Actual (total) program run time: " << (double)(end_t-start_t)/CLOCKS_PER_SEC << "\n";
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
  if (aveKE) {
    cout << "KE=" << average(simulation.getStatistic(0)) << ";\n";
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
    cout << "Print[\"Average Profile\"]\nListLinePlot[aveProf,PlotStyle->Black,ImageSize->Large]\n";
  }

  return 0;
}
