/* 
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#include "../control/Creator.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"

#include "ArgParse.h"

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

  // Set up a data recorder
  DataRecord *dataRecord = new DataRecord;
  integrator.setDataRecord(dataRecord);
  dataRecord->addStatFunction(StatFunc_AveKE, "KE");

  // Set parameters
  double time = 60.;

  // Can get the time from the command line
  ArgParse parser;
  try {
    parser.set(argc, argv);
  }
  catch (ArgParse::IllegalToken token) {
    cout << "Illegal Token: " << token.c << ". Exiting.\n";
    exit(1);
  }
  parser.get("time", time);
  // Make sure we didn't enter any illegal tokens (ones not listed above) on the command line
   try {
     parser.check();
   }
   catch (ArgParse::UncheckedToken illegal) {
     cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
     exit(1);
   }

  integrator.initialize(time);

  // Print initial message
  cout << "Starting integration.\n";

  // Run here
  integrator.integrate();

  // Print a final message
  cout << "Integration ended.\n";
  double runTime = dataRecord->getElapsedTime();
  cout << "Run time: " << runTime << endl;
  cout << "Ratio: " << time / runTime << endl;

  // Get data from the data recorder
  dataRecord->writeData("RunData", simData);
  
  int numStatFuncs = dataRecord->getNumberOfStatFunctions();
  for (int i=0; i<numStatFuncs; ++i) {
    auto data   = dataRecord->getStatFunctionData(i);
    string name = dataRecord->getStatFunctionName(i);
    cout << name << "=" << data << endl;
  }
      
  
#ifdef USE_MPI
  MPI::Finalize();
#endif

  // Clean up
  delete simData;
  delete dataRecord;

  return 0;
}
