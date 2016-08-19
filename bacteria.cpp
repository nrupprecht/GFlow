#include "Simulator.h"

int main(int argc, char** argv) {
  auto start_t = clock();
  // Parameters
  double width = 2.;
  double height = 1.;
  double radius = 0.02;
  double velocity = 0.;
  double temperature = 8.;
  double phi = 0.001;
  double replenish = 0.;
  double time = 600.;
  double start = 0;

  // Display parameters
  bool animate = true;
  bool dispProfile = false;
  bool dispAveProfile = false;
  bool recFields = true;

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
  opt = parser.find("temperature");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> temperature;
  }
  opt = parser.find("phi");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> phi;
  }
  opt = parser.find("replenish");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> replenish;
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
  opt = parser.find("animate");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> animate;
  }
  opt = parser.find("recFields");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> recFields;
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
  //----------------------------------------

  // Dependent variables
  double Vol = width*height;
  int number = Vol/(PI*sqr(radius))*phi;
  
  // Seed random number generators
  srand48( std::time(0) );
  srand( std::time(0) );
  
  //----------------------------------------
  
  Simulator simulation;
  simulation.setStartRecording(start);
  simulation.createBacteriaBox(number, radius, width, height, velocity);
  simulation.setReplenish(replenish);
  simulation.setTemperature(temperature);
  simulation.setRecFields(recFields);
  simulation.bacteriaRun(time);
  auto end_t = clock();
  
  /// Print condition summary
  cout << "Dimensions: " << width << " x " << height << "\n";
  cout << "Radius: " << radius << "\n";
  cout << "Fluid Velocity: " << velocity << "\n";
  cout << "Phi: " << phi << ", Number: " << number << ", (Actual Phi: " << number*PI*sqr(radius)/(width*height) << ")\n";
  cout << "Sim Time: " << time << ", Run time: " << simulation.getRunTime() << ", Ratio: " << time/simulation.getRunTime() << endl;
  cout << "Start Time: " << start << "\n";
  cout << "Actual (total) program run time: " << (double)(end_t-start_t)/CLOCKS_PER_SEC << "\n";
  cout << "Iters: " << simulation.getIter() << "\n\n";
  
  /// Print data
  if (animate) {
    cout << mmPreproc(simulation.printWatchList()) << endl;
    cout << "walls=" << simulation.printWalls() << ";\n";
    cout << simulation.printAnimationCommand() << endl;
  }
  
  if (dispProfile) {
    cout << "prof=" << simulation.getProfile() << ";\n";
    cout << "Print[\"Profile\"]\nMatrixPlot[prof,ImageSize->Large]\n"; 
  }
  if (dispAveProfile) {
    cout << "aveProf=" << simulation.getAveProfile() << ";\n";
    cout << "Print[\"Average Profile\"]\nListLinePlot[aveProf,PlotStyle->Black,ImageSize->Large]\n";
  }

  // Print Fields //**
  if (recFields) {
    cout << "res=" << mmPreproc(simulation.printResourceRec()) << ";\n";
    cout << "Print[\"Resource\"]\n";
    cout << "resframes=Table[MatrixPlot[res[[i]]],{i,1,Length[res]}];\n";
    cout << "ListAnimate[resframes]\n";
    cout << "wst=" << mmPreproc(simulation.printWasteRec()) << ";\n";
    cout << "Print[\"Waste\"]\n";
    cout << "wstframes=Table[MatrixPlot[wst[[i]]],{i,1,Length[wst]}];\n";
    cout << "ListAnimate[wstframes]\n";
    cout << "fit=" << mmPreproc(simulation.printFitnessRec()) << ";\n";
    cout << "Print[\"Fitness\"]\n";
    cout << "fitframes=Table[MatrixPlot[fit[[i]]],{i,1,Length[fit]}];\n";
    cout << "ListAnimate[fitframes]\n";  
  }

  return 0;
}
