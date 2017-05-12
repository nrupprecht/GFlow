/* 
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

// Includes
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "../integrators/VelocityVerletIntegrator.hpp"

using namespace GFlow;

int main (int argc, char** argv) {

  int rank = 0, numProc = 0;
  // Set up MPI
#ifdef USE_MPI
  MPI::Init();
  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();
#endif

  // Create from file or command line args
  VelocityVerletIntegrator integrator;

  // Run here
  integrator.integrate();

  // Print a final message
  cout << "Integration ended.\n";
  
#ifdef USE_MPI
  MPI::Finalize();
#endif

  return 0;
}
