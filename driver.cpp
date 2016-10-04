#include "Simulator.h"

int main(int argc, char** argv) {
  auto start_t = clock();
  // Parameters
  double width = 4.;     // Height of the pipe
  double height = 2.;    // Width of the pipe
  double radius = 0.05;  // Disc radius
  double velocity = 0.5; // Fluid velocity (at center)
  double phi = 0.5;      // Packing density
  double time = 600.;    // How long the entire simulation should last
  double start = 30;     // At what time we start recording data
  double pA = 0.;        // What percent of the particles are active
  double activeF = 0.25; // Active force (default is 5)
  int samplePoints = -1; // How many bins we should use when making a density profile

  // Display parameters
  bool animate = false;
  bool dispKE = false;
  bool aveKE = false;
  bool dispFlow = false;
  bool dispProfile = false;
  bool dispAveProfile = true;
  bool dispVelDist = true;

  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser(argc, argv);
  pair<string,string> opt;
  stringstream stream;
  parser.get("width", width);
  parser.get("height", height);
  parser.get("radius", radius);
  parser.get("velocity", velocity);
  parser.get("phi", phi);
  parser.get("time", time);
  parser.get("start", start);
  parser.get("active", pA);
  parser.get("force", activeF);
  parser.get("points", samplePoints);
  parser.get("animate", animate);
  parser.get("dispKE", dispKE);
  parser.get("aveKE", aveKE);
  parser.get("flow", dispFlow);
  parser.get("profileMap", dispProfile);
  parser.get("profile", dispAveProfile);
  parser.get("velDist", dispVelDist);
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
  cout << "Command: ";
  for (int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << "\n-------------------------------------\n";

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
  if (dispVelDist) {
    cout << "velDist=" << simulation.getVelocityDistribution() << ";\n";
    cout << "Print[\"Velocity Distribution\"]\nListLinePlot[velDist,PlotStyle->Black,ImageSize->Large,PlotRange->All]\n";
    cout << "aVelDist=" << simulation.getAuxVelocityDistribution() << ";\n";
    cout << "Print[\"Auxilary Velocity Distribution\"]\nListLinePlot[aVelDist,PlotStyle->Black,ImageSize->Large,PlotRange->All]\n";
  }

  return 0;
}
