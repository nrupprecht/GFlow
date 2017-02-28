#include "Simulator.h"

int main(int argc, char* argv[]) {
  // Parameters
  double gap = 0.2; // Hopper gap
  double width = 3;  // Width of the hopper
  double height = 4; // Height of the hopper
  int number = 500;  // Number of particles in the hopper
  int timeIters = 25;    // The number of time points to use
  int trials = 10;   // Number of trials to do
  double radius = 0.05; // Radius of the particles
  double actPer = 0.; // Percent of the particles that are active
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
  vector<vect<> > data; // use x = time, y = % jammed
  double minTime = 0, maxTime = 10;
  double dt = (maxTime-minTime)/(timeIters+1);
  for (int i=0; i<timeIters; i++) {
    int fullTime = 0; // Number of simulations that didn't jam in the given time
    double time = (i+1)*dt+minTime;
    for (int i=0; i<trials; i++) {
      simulation.createHopper(number, radius, gap, width, height, actPer);
      simulation.run(time+delayTime+3.0); // 3 - when the temp wall is taken away, delayTime - how long it takes to detect a jam
      if (!simulation.getDelayTriggeredExit()) fullTime++;
    }
    data.push_back(vect<>(time,(double)fullTime/trials));
  }
  clock_t end = clock();

  // Print data
  cout << endl;
  cout << data << endl;
  cout << endl;
  cout << "Number: " << number << ", Trials: " << trials << endl;
  cout << "Active: " << actPer << endl;
  cout << "Gap: " << gap << ", Radius: " << radius << ", Max time: " << time << endl;
  cout << "Run time: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
  cout << endl;

  return 0;
}
