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
  double radius = 0.05;
  double time = 1.;

  // Initialize MPI
  int rank, numProc;
  MPI::Init(argc, argv);
  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();
  if (rank==0) cout << "Using " << numProc << " processors.\n";
  MPI::COMM_WORLD.Barrier();

  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser(argc, argv);
  parser.get("number", number);
  parser.get("radius", radius);
  parser.get("time", time);
  //----------------------------------------

  // Set up the simulation
  GFlowBase simulator;  
  simulator.createSquare(number, radius);
  simulator.run(time);

  /// Print condition summary
  if (rank==0) {
    cout << "----------------------- RUN SUMMARY -----------------------\n";
    cout << "Command: ";
    for (int i=0; i<argc; i++) cout << argv[i] << " ";
    cout << endl << endl; // Line break
    cout << "Number: " << number << endl;
    cout << "Sim Time: " << time << ", Run Time: " << simulator.getRunTime() << endl;
    cout << "\n----------------------- END SUMMARY -----------------------\n\n"; 

    /// Print recorded data
    cout << "pos=" << simulator.getPositionRecord() << ";\n";

  }
  // End MPI
  MPI::Finalize();
  return 0;
}
