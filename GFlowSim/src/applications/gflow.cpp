/* 
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#include "../creation/Creator.hpp"
#include "../creation/FileParser.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../../include/ArgParse.h"

using namespace GFlow;

int main (int argc, char** argv) {
  
  string config = "";
  string writeDirectory = "";
  string loadFile = "";
  string saveFile = "";
  RealType time = 10;
  // Record options
  bool animate = false;
  int recOption = 0;
  bool recPerf = false;
  bool recMvRatio = false;
  bool recDt = false;
  bool recDelay = false;
  int fps = 15;
  // Stat options
  bool maxRatio = false;
  bool ratio = false;
  bool KE = false;
  bool maxV = false;
  bool aveV = false;
  bool aveF = false;
  bool maxF = false;
  
  bool print = false;       // Whether we should print stat data to the screen
  bool quiet = false;
  // Performance options
  int adjust = -1;       // Whether to auto-adjust the time step
  int adjustDelay = -1; // Whether to auto-adjust the update delay
  RealType maxDt = -1;
  RealType minDt = -1;
  RealType dt = -1;

  // Can get options from the command line
  ArgParse parser;
  try {
    parser.set(argc, argv);
  }
  catch (ArgParse::IllegalToken token) {
    cout << "Illegal Argument: " << token.c << ". Exiting.\n";
    exit(1);
  }
  // Logistics options
  parser.get("config", config);
  parser.get("writeDirectory", writeDirectory); 
  parser.get("loadFile", loadFile); 
  parser.get("saveFile", saveFile);
  parser.get("time", time);
  // Animation options
  parser.get("animate", animate);
  parser.get("recOption", recOption);
  parser.get("recPerf", recPerf);
  parser.get("recMvRatio", recMvRatio);
  parser.get("recDt", recDt);
  parser.get("recDelay", recDelay);
  parser.get("fps", fps);
  parser.get("maxRatio", maxRatio);
  parser.get("ratio", ratio);
  parser.get("KE", KE);
  parser.get("maxV", maxV);
  parser.get("aveV", aveV);
  parser.get("maxF", maxF);
  parser.get("aveF", aveF);
  parser.get("print", print);
  parser.get("quiet", quiet);
  // Performance options
  parser.get("adjust", adjust);
  parser.get("adjustDelay", adjustDelay);
  parser.get("maxDt", maxDt);
  parser.get("minDt", minDt);
  parser.get("dt", dt);
  // Make sure we didn't enter any illegal tokens (ones not listed above) on the command line
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }
  
  // Set up MPI
#if USE_MPI == 1
#if _CLANG_ == 1
  int rank, numProc;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
#else
  MPI::Init(argc, argv);
  int rank = MPI::COMM_WORLD.Get_rank();
  int numProc = MPI::COMM_WORLD.Get_size();
#endif
#endif

  // Create a data record
  DataRecord *dataRecord = new DataRecord;
  
  // Integrator
  Integrator *integrator = nullptr;
  
  // Create from file or command line args
  SimData *simData = nullptr;
  if (config!="") {
    FileParser fileParser;
    fileParser.setDataRecord(dataRecord);
    try {
      fileParser.parse(config, simData, integrator);
    }
    catch (FileParser::FileDoesNotExist file) {
      cout << "File [" << file.name << "] does not exist. Trying [samples/" << file.name << "].\n";
      try {
	fileParser.parse("samples/"+config, simData, integrator);
      }
      catch (FileParser::FileDoesNotExist file) {
	cout << "File [" << file.name << "] also does not exist. Exiting.\n";
	exit(0);
      }
      cout << "You're lucky, [samples/" << file.name << "] does exist. Assuming you meant that file.\n";
    }
    catch (FileParser::UnrecognizedToken token) {
      cout << "Unrecognized token [" << token.token << "]. Exiting.\n";
      exit(0);
    }
    if (simData==nullptr) {
      cout << "Error occured while parsing, sim data is null. Exiting.\n";
      exit(0);
    }
  }
  else if (loadFile!="") {
    FileParser fileParser;
    fileParser.setDataRecord(dataRecord);
    try {
      simData = fileParser.loadFromFile(loadFile);
    }
    catch (FileParser::FileDoesNotExist file) {
      cout << "File [" << file.name << "] does not exist. Trying [samples/" << file.name << "]. Exiting.\n";
      exit(0);
    }
    catch (GFlow::BadIstreamRead read) {
      cout << "Parsing failed when it attempted to read in a [" << read.type << "], expected [" << printChar(read.expected) << "], found [" << printChar(read.unexpected) << "]. Exiting.\n";
      exit(0);
    }
  }
  else {
    Creator creator;
    simData = creator.create();
  }

  // Make sure simData is non-null
  if (simData==nullptr) {
    cout << "SimData is null. Exiting." << endl;
    exit(0);
  }

  // Create an integrator
  // VelocityVerletIntegrator integrator(simData);
  if (maxDt>0) reinterpret_cast<VelocityVerletIntegrator*>(integrator)->setMaxTimeStep(maxDt);
  if (minDt>0) reinterpret_cast<VelocityVerletIntegrator*>(integrator)->setMinTimeStep(minDt);
  if (dt>0)    integrator->setDt(dt);

  // Set up the data recorder
  if (writeDirectory!="") dataRecord->setWriteDirectory(writeDirectory);
  integrator->setDataRecord(dataRecord);
  // Set record options
  if (dataRecord) {
    // Set statistic options
    if (KE)       dataRecord->addStatFunction(StatFunc_AveKE, "KE");
    if (maxRatio) dataRecord->addStatFunction(StatFunc_MaxVelocitySigmaRatio, "maxRatio");
    if (ratio)    dataRecord->addStatFunction(StatFunc_MinSigmaVelocityRatio, "ratio");
    if (maxV)     dataRecord->addStatFunction(StatFunc_MaxSpeed, "maxV");
    if (aveV)     dataRecord->addStatFunction(StatFunc_AveSpeed, "aveV");
    if (maxF)     dataRecord->addStatFunction(StatFunc_MaxForce, "maxF");
    if (aveF)     dataRecord->addStatFunction(StatFunc_AveForce, "aveF");
    // Set recording options
    dataRecord->setRecPos(animate);
    dataRecord->setRecOption(recOption);
    dataRecord->setRecPerf(recPerf);
    dataRecord->setRecMvRatio(recMvRatio);
    dataRecord->setRecDt(recDt);
    dataRecord->setRecDelay(recDelay);
    dataRecord->setDelay(1./static_cast<RealType>(fps));
  }

  // Set run time
  integrator->initialize(time);
  if (adjust>-1) reinterpret_cast<VelocityVerletIntegrator*>(integrator)->setAdjustTimeStep(adjust);
  if (adjustDelay>-1) reinterpret_cast<VelocityVerletIntegrator*>(integrator)->setAdjustUpdateDelay(adjustDelay);

  // Print initial message
#if USE_MPI == 1
  if (rank==0) {
#endif
    if (!quiet) cout << "Starting integration.\n";
#if USE_MPI == 1
  }
#endif

  // Run the integrator
  integrator->integrate();

  // Print a final message
  #if USE_MPI == 1
  if (rank==0) {
#endif
    if (!quiet) cout << "Integration ended.\n";
#if USE_MPI == 1
  }
#endif

#if USE_MPI == 1
  if (rank==0) {
#endif
    if (dataRecord) {
      // Print out time and ratio
      double runTime = dataRecord->getElapsedTime();
      if (!quiet) {
	cout << "Run time: " << runTime << endl;
	cout << "Ratio: " << time / runTime << endl;
      }
      
      // Write animation data to files
      dataRecord->writeData(simData);
      
      // Write sectorization data to file
      dataRecord->writeRunSummary(simData, integrator);
      
      // Write out stat function data - for now
      if (print) {
	int numStatFuncs = dataRecord->getNumberOfStatFunctions();
	for (int i=0; i<numStatFuncs; ++i) {
	  auto data   = dataRecord->getStatFunctionData(i);
	  string name = dataRecord->getStatFunctionName(i);
	  cout << name << "=" << mmPreproc(data) << ";\n";
	  cout << "ListLinePlot[" << name << ",ImageSize->Large,PlotStyle->Black]\n";
	}
      }
    }
#if USE_MPI == 1
  }
#endif

  // Save arrangement if requested
#if USE_MPI == 1
  if (rank==0) {
#endif
    if (saveFile!="") {
      FileParser fileParser;
      fileParser.saveToFile(simData, saveFile);
    }
#if USE_MPI == 1
  }
#endif
  
  // Finalize mpi
#if USE_MPI == 1
#if _CLANG_ == 1
  MPI_Finalize();
#else
  MPI::Finalize();
#endif
#endif

  // Clean up
  if (simData)    delete simData;
  if (dataRecord) delete dataRecord;

  // The program is done. Return.
  return 0;
}
