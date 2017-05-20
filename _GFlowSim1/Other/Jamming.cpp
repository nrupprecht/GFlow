#include "Simulator.h"

/// Simulates Hopper jamming

int main(int argc, char* argv[]) {
  // Parameters
  double gap = 0.15; // Hopper gap
  double width = 3;  // Width of the hopper
  double height = 4; // Height of the hopper
  int number = 500;  // Number of particles in the hopper
  int trials = 100;   // Number of trials to do
  double time = 180; // Max time we let the simulation run for
  double radius = 0.05; // Radius of the particles
  double actPer = 0.1; // Percent of the particles that are active

  double delayTime = 10.; // How long to wait before we say the system is jammed

  // First parameter: gap width (optional)
  stringstream stream;
  if (argc>1) {
    stream << argv[1];
    stream >> gap;
  }

  srand48(std::time(0)); // Seed for drand48
  srand(std::time(0));   // Seed for rand (I don't think this is neccessary)

  Simulator simulation;

  simulation.setMarkWatch(true);
  simulation.setDelayTime(delayTime);
  simulation.setStartTime(3.);

  clock_t start = clock();
  int fullTime = 0; // Number of simulations that didn't jam in the given time
  vector<double> jamTimes, watchSlopes;
  vector<int> markSize, size;
  for (int i=0; i<trials; i++) {
    simulation.createHopper(number, radius, gap, width, height, actPer);
    simulation.run(time);    
    double T = simulation.getTime();
    size.push_back(simulation.getPSize());
    if (!simulation.getDelayTriggeredExit()) fullTime++;
    double jtime = simulation.getMarkDiff();
    jamTimes.push_back(jtime);
    watchSlopes.push_back(simulation.getMarkSlope());
    markSize.push_back(simulation.getMarkSize());
  }
  clock_t end = clock();

  // Find the average
  double aveTime = 0, aveSlope, aveSize;
  for (auto t : jamTimes) aveTime += t;
  aveTime /= trials;
  for (auto s : watchSlopes) aveSlope += s;
  aveSlope /= trials;
  for (auto s : markSize) aveSize += s;
  aveSize /= trials;

  // Find the standard deviations
  double stdDevA = 0, stdDevB = 0, stdDevC = 0;
  for (auto t : jamTimes) stdDevA += sqr(aveTime-t);
  stdDevA /= trials;
  stdDevA = sqrt(stdDevA);
  for (auto t : watchSlopes) stdDevB += sqr(aveSlope-t);
  stdDevB /= trials;
  stdDevB = sqrt(stdDevB);
  for (auto t : markSize) stdDevC += sqr(aveSize-t);
  stdDevC /= trials;
  stdDevC = sqrt(stdDevC);

  // Print data
  cout << endl;
  cout << "jamTimes=" << jamTimes << ";\n";
  cout << "sizes=" << size << ";\n";
  cout << endl;
  cout << "Number: " << number << ", Trials: " << trials << endl;
  cout << "Active: " << actPer << endl;
  cout << "Gap: " << gap << ", Radius: " << radius << ", Max time: " << time << endl;
  cout << "Didn't jam " << fullTime << " times.\n";
  cout << "Ave jam time: " << aveTime << ", std dev: " << stdDevA << endl;
  cout << "Slope: " << aveSlope << ", std dev: " << stdDevB << endl;
  cout << "Number of marks: " << aveSize << ", std dev: " << stdDevB << endl;
  cout << "Run time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
  cout << endl;

  return 0;
}
