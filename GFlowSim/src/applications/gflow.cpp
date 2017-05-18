/* 
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#include "../creation/Creator.hpp"
#include "../creation/FileParser.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"

#include "ArgParse.h"

using namespace GFlow;

int main (int argc, char** argv) {

  string config = "";
  RealType time = 10;
  bool animate = false;
  bool ratio = false;
  bool KE = false;
  bool perf = false;
  bool adjust = true;       // Whether to auto-adjust the time step
  bool updateAdjust = true; // Whether to auto-adjust the update delay

  // Can get options from the command line
  ArgParse parser;
  try {
    parser.set(argc, argv);
  }
  catch (ArgParse::IllegalToken token) {
    cout << "Illegal Token: " << token.c << ". Exiting.\n";
    exit(1);
  }
  parser.get("config", config);
  parser.get("time", time);
  parser.get("animate", animate);
  parser.get("ratio", ratio);
  parser.get("KE", KE);
  parser.get("perf", perf);
  parser.get("adjust", adjust);
  parser.get("updateAdjust", updateAdjust);
  // Make sure we didn't enter any illegal tokens (ones not listed above) on the command line
   try {
     parser.check();
   }
   catch (ArgParse::UncheckedToken illegal) {
     cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
     exit(1);
   }

  // Set up MPI
#ifdef USE_MPI
  MPI::Init();
  int rank = MPI::COMM_WORLD.Get_rank();
  int numProc = MPI::COMM_WORLD.Get_size();
#endif

  // Create from file or command line args
  Creator simCreator;
  SimData *simData = nullptr;
  if (config!="") {
    FileParser fileParser;
    try {
      simData = fileParser.parse(config);  
    }
    catch (FileParser::FileDoesNotExist file) {
      cout << "File [" << file.name << "] does not exist. Exiting.\n";
      exit(0);
    }
    if (simData==nullptr) {
      cout << "Error occured while parsing. Exiting.\n";
      exit(0);
    }
  }
  else simData = simCreator.create();

  // Make sure simData is non-null
  if (simData==nullptr) {
    cout << "SimData is null. Exiting." << endl;
    exit(0);
  }

  // Create an integrator
  VelocityVerletIntegrator integrator(simData);

  // Set up a data recorder
  DataRecord *dataRecord = new DataRecord;
  integrator.setDataRecord(dataRecord);
  // Set record options
  if (dataRecord) {
    if (KE) dataRecord->addStatFunction(StatFunc_AveKE, "KE");
    if (ratio) dataRecord->addStatFunction(StatFunc_MaxVelocitySigmaRatio, "ratio");

    dataRecord->setRecPos(animate);
    dataRecord->setRecPerf(perf);
  }

  // Set run time
  integrator.initialize(time);
  integrator.setAdjustTimeStep(adjust);
  integrator.setAdjustUpdateDelay(updateAdjust);

  // Print initial message
  cout << "Starting integration.\n";

  // Run the integrator
  integrator.integrate();

  // Print a final message
  cout << "Integration ended.\n";
  if (dataRecord) {
    // Print out time and ratio
    double runTime = dataRecord->getElapsedTime();
    cout << "Run time: " << runTime << endl;
    cout << "Ratio: " << time / runTime << endl;

    // Write animation data to files
    dataRecord->writeData(simData);
    
    // Write sectorization data to file
    dataRecord->writeRunSummary(simData, &integrator);

    // Write out stat function data - for now
    int numStatFuncs = dataRecord->getNumberOfStatFunctions();
    for (int i=0; i<numStatFuncs; ++i) {
      auto data   = dataRecord->getStatFunctionData(i);
      string name = dataRecord->getStatFunctionName(i);
      cout << name << "=" << mmPreproc(data,3) << ";\n";
      cout << "ListLinePlot[" << name << ",ImageSize->Large,PlotStyle->Black]\n";
    }

    // Print performance
    auto performanceRecord = dataRecord->getPerformanceRecord();
    if (performanceRecord.size()>0) {
      cout << "perf=" << mmPreproc(performanceRecord,3) << ";\n";
      cout << "ListLinePlot[perf,ImageSize->Large,PlotStyle->Black]\n";
    }
  }
      
  
#ifdef USE_MPI
  MPI::Finalize();
#endif

  // Clean up
  if (simData)    delete simData;
  if (dataRecord) delete dataRecord;

  return 0;
}
