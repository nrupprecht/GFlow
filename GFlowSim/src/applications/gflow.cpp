/* 
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#include "../control/Creator.hpp"
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
  Creator simCreator;
  SimData *simData = simCreator.create();
  VelocityVerletIntegrator integrator(simData);
  DataRecord *dataRecord = new DataRecord;
  integrator.setDataRecord(dataRecord);
  
  // Set parameters
  integrator.initialize(5.);

  // Print initial message
  cout << "Starting integration.\n";

  // Run here
  integrator.integrate();

  // Print a final message
  cout << "Integration ended.\n";

  // Write data
  dataRecord->writeData("RunData", simData);
  
#ifdef USE_MPI
  MPI::Finalize();
#endif

  // Clean up
  delete simData;
  delete dataRecord;

  return 0;
}
