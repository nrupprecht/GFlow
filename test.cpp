#include "TestSec.h"
#include "Simulator.h"

int main(int argv, char** argc) {
  // Set up
  TestSec sectorization(0,1,0,1);
  double epsilon = 1e-4;
  int N = 1000;
  double radius = 0.01;
  sectorization.setSectorDim(2*radius);
  // sectorization.setParticlesInteract(false); //** For now
  // For display purposes
  vector<vector<vect<> > > positions;
  // Add some particles
  for (int i=0; i<N; ++i) {
    // Random position
    vect<> pos(drand48(), drand48());
    // Random velocity
    double angle = 2*PI*drand48();
    vect<> velocity = 0.1*drand48()*vect<>(cos(angle), sin(angle));
    sectorization.addParticle(pos, radius, velocity);
  }
  
  // Run for a while
  double runLength = 30.;
  double time = 0, dispDelay = 1./15, lastDisp = -dispDelay;
  double rearrangeDelay = 0.5, lastRearrange = 0;
  int iter = 0;
  sectorization.rearrangeParticleList();
  // Start the clock
  auto start_t = clock();
  while (time<runLength) {
    sectorization.interactions();
    sectorization.updateParticles(epsilon);
    sectorization.updateSectors();
    // Record display data
    /*
     if (time-lastDisp>=dispDelay) {
       positions.push_back(sectorization.getPositions());
       lastDisp = time;
     }
    */
    if (time-lastRearrange>rearrangeDelay) {
      sectorization.rearrangeParticleList();
      lastRearrange = time;
    }
    
    time += epsilon;
    iter++;
  }
  
  auto end_t = clock(); // End timing
  double realTime = (double)(end_t-start_t)/CLOCKS_PER_SEC;
  cout << "Real Time: " << realTime << ", Iterations: " << iter << endl;

  /*
  Simulator simulator;
  simulator.createSquare(N, 0, radius, 1, 1);
  simulator.setParticleInteraction(false);
  simulator.run(runLength);
  cout << "Simulator Real Time: " << simulator.getRunTime() << endl;
  */
  
  /*
  cout << "pos=" << mmPreproc(positions) << ";\n";
  cout << "R0=" << radius << ";\nlen=Length[pos];\nscale=100;\nG0=Table[Graphics[Table[{Black,Circle[pos[[i]][[j]],R0]},{j,1,Length[pos[[i]]]}]],{i,1,len}];\nframes=Table[Show[G0[[i]],PlotRange->{{0,1},{0,1}},ImageSize->{scale*4,scale*4}],{i,1,len}];\nvid=ListAnimate[frames,AnimationRate->15]";
  */
  return 0;
}
