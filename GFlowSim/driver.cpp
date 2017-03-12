/// driver.cpp --- Driver file
/// Nathaniel Rupprecht 2016
///

#include "../ArgParse.h"
#include "GFlowBase.h"
//#include "GFlow.h"

using MPI::COMM_WORLD;

#include <iostream>
using std::cout;
using std::endl;

int main(int argc, char** argv) {
  // Simulation parameters
  int number = 100;
  double width = 4;
  double height = 4;
  double radius = 0.05;
  double bR = 0.2;
  double density = 10;
  double velocity = 0.25;
  double dispersion = 0;
  double time = 1.;
  double phi = -1;
  bool interact = true;

  // Animation Paramaters
  bool animate = false;
  bool KE      = false;
  bool novid   = false;

  // Simulation type
  bool square = true;
  bool buoyancy = false;

  // Initialize MPI
  int rank, numProc;
  MPI::Init();
  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();
  MPI::COMM_WORLD.Barrier();

  // Seed random number generators
  srand48( std::time(0) );
  srand( std::time(0) );
  seedNormalDistribution();

  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser(argc, argv);
  parser.get("number", number);
  parser.get("width", width);
  parser.get("height", height);
  parser.get("radius", radius);
  parser.get("bR", bR);
  parser.get("density", density);
  parser.get("velocity", velocity);
  parser.get("dispersion", dispersion);
  parser.get("time", time);
  parser.get("phi", phi);
  parser.get("interact", interact);
  parser.get("animate", animate);
  parser.get("KE", KE);
  parser.get("novid", novid);
  parser.get("square", square);
  parser.get("buoyancy", buoyancy);
  //----------------------------------------

  // Calculate number of particles given a packing fraction
  if (phi!=-1) {
    double Vol = width*height;
    number = Vol/(PI*sqr(radius))*phi;
  }

  /// Print condition summary
  if (rank==0) {
    cout << "----------------------- RUN SUMMARY -----------------------\n";
    cout << "Command: ";
    for (int i=0; i<argc; i++) cout << argv[i] << " ";
    cout << endl << endl; // Line break
    cout << "Using " << numProc << " processors.\n";
    cout << "Number: " << number << endl;
    cout << "  ..........................\n";
  }
  MPI::COMM_WORLD.Barrier();

  // Set up the simulation
  GFlowBase simulator;  
  if (buoyancy) simulator.createBuoyancyBox(radius, bR, density, width, height, velocity, dispersion);
  else if (square) simulator.createSquare(number, radius, width, height, velocity, dispersion);
  else throw false; // No selection
  simulator.setDoInteractions(interact);
  simulator.setRecPositions(animate);
  simulator.setRecKE(KE);

  if (rank==0) {
    cout << "Dimensions: " << simulator.getWidth() << " x " << simulator.getHeight() << endl;
    cout << "Set up time: " << simulator.getSetUpTime() << endl;
    cout << "  ..........................\n";
  }

  // Run the simulation
  simulator.run(time);

  // Print Run summary
  if (rank==0) {
    double runTime = simulator.getRunTime(), transferTime = simulator.getTransferTime();
    int iters = simulator.getIter(), ndx = simulator.getNDX(), ndy = simulator.getNDY();
    cout << "Domains: " << ndx << " x " << ndy << ", Total: " << ndx*ndy << endl;
    cout << "Run Time: " << runTime << ", Sim Time: " << time << endl;
    cout << "Iterations: " << iters << endl;
    cout << "Transfer Time: " << transferTime << " (" << (runTime>0 ? toStr(transferTime/runTime*100) : "---") << "%)" << endl;
    cout << "Ratio: " << (runTime>0 ? toStr(time/runTime) : "---") << "\tTime per iter per particle: " << (iters>0 && number>0 ? mmPreproc(time/iters/number) : "---");
    cout << "\n----------------------- END SUMMARY -----------------------\n\n"; 

    /// Print recorded data
    if (animate) cout << simulator.printAnimationCommand(novid) << endl;
    if (KE) {
      cout << "ke=" << mmPreproc(simulator.getKERecord()) << ";\n";
      cout << "ListLinePlot[ke,ImageSize->Large,PlotStyle->Black,PlotRange->{0,1.1*Max[ke]}]\n";
    }

  }
  // End MPI
  MPI::Finalize();
  return 0;
}
