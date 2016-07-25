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
  double width = 5.;
  double height = 2;
  int number = 250; // Total number of particles
  double pA = 0.1;  // Portion that are active
  int trials = 1;   // Number of trials to average over
  double time = 30; // Time to run the simulation for
  double startRec = 3; // What time to start recording data
  double radius = 0.05;
  double rA = 0.05;
  
  double nA = pA*number;
  
  // Seed random number generators
  srand48( std::time(0) );
  srand( std::time(0) );

  //----------------------------------------

  Simulator simulation;
  
  simulation.addStatistic(statPassiveFlow);
  simulation.addStatistic(statActiveFlow);

  simulation.setStartRecording(startRec);
  auto start = clock();
  vector<vect<> > dataA, dataP;
  for (int i=0; i<6; i++) {
    double B = 0.02*i; 
    vect<> bias(B, 0);
    double act=0, pass=0;
    for (int i=0; i<trials; i++) {
      simulation.createControlPipe(number-nA, nA, radius, 0, bias, rA, width, height);
      pass += average(simulation.getStatistic(0));
      act += average(simulation.getStatistic(1));
      simulation.run(time);
    }
    act /= (double)trials;
    pass /= (double)trials;
    
    dataA.push_back(vect<>(B, act));
    dataP.push_back(vect<>(B, pass));
  }
  auto end = clock();
  
  int psize = simulation.getPSize(), asize = simulation.getASize();
  cout << "Passive: " << psize << ", Active: " << asize << endl;
  cout << "Pipe width, height: " << width << ", " << height << endl;
  cout << "Packing: " << (psize*PI*sqr(radius)+asize*PI*sqr(rA))/(width*height) << endl;
  cout << "Time: " << (double)(end-start)/CLOCKS_PER_SEC;

  cout << "act=" << dataA << ";\n";
  cout << "pass=" << dataP << ";";
  
  return 0;
}
