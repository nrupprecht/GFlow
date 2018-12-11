// For argument parsing
#include "../utility/ArgParse.hpp"

// Simulation creators
#include "../allcreators.hpp"

// All data objects to choose from
#include "../alldataobjects.hpp"

// All the modifiers
#include "../allmodifiers.hpp"

/*
*  --- NOTES:
*  Running functions dynamically:
*  <https://stackoverflow.com/questions/11016078/is-it-possible-to-create-a-function-dynamically-during-runtime-in-c>
*  <http://burnttoys.blogspot.com/2011/04/how-to-allocate-executable-memory-on.html>
*
*/

using namespace GFlowSimulation;

int main(int argc, char **argv) {
  // --- MPI
  int rank(0), numProc(1);

  #if USE_MPI == 1
  #if _CLANG_ == 1
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  #else
  MPI::Init(argc, argv);
  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();
  #endif
  cout << "Initialized MPI. Rank " << rank << "\n";
  #endif

  // --- Options

  // Type of simulation
  bool bipartite_flag = false;
  bool debug_flag = false;
  string load = "";

  // Data to gather
  bool animate = false; // Record positions
  bool snapshot = false;
  bool sectorData = false; // Record sector information
  bool ke = false; // Record kinetic energy
  bool keTypes = false;
  bool totalKE = false; // Record average kinetic energy (per particle)
  bool secRemake = false; 
  bool bdForces = false;
  bool timestep = false;
  bool averages = false;
  bool aveV = false;
  bool aveP = false;
  bool minDistances = false;
  bool percolation = false;
  bool psnapshot = false;
  bool memdist = false;
  bool pressure = false;
  bool numberdata = false;
  
  // Other options
  int dimensions = 2;
  RealType skin = 0.;
  bool quiet = false;
  RealType gravity = 0.;
  bool damping = false;
  bool adjustDT = false;
  int target_steps = -1;
  int step_delay = -1;
  RealType startRecTime = 0;
  RealType fps = -1.;
  RealType videoLength = -1.;
  RealType dt = 0.0001;
  RealType time = 10.;
  string writeDirectory = "RunData";
  int boundary = 1;
  string monitor = ""; // Monitor file

  // Modifiers
  RealType temperature = 0;

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("bipartite", bipartite_flag);
  parser.get("debug", debug_flag); 
  parser.get("load", load);
  parser.get("animate", animate);
  parser.get("snapshot", snapshot);
  parser.get("sectorData", sectorData);
  parser.get("KE", ke);
  parser.get("KETypes", keTypes);
  parser.get("totalKE", totalKE);
  parser.get("secRemake", secRemake);
  parser.get("bdForces", bdForces);
  parser.get("timestep", timestep);
  parser.get("averages", averages);
  parser.get("aveV", aveV);
  parser.get("aveP", aveP);
  parser.get("minDistances", minDistances);
  parser.get("percolation", percolation);
  parser.get("psnapshot", psnapshot);
  parser.get("memdist", memdist);
  parser.get("pressure", pressure);
  parser.get("numberdata", numberdata);
  parser.get("dimensions", dimensions);
  parser.get("skin", skin);
  parser.get("quiet", quiet);
  parser.get("gravity", gravity);
  parser.get("damping", damping);
  parser.get("adjustDT", adjustDT);
  parser.get("target_steps", target_steps);
  parser.get("step_delay", step_delay);
  parser.get("startRec", startRecTime);
  parser.get("fps", fps);
  parser.get("videoLength", videoLength);
  parser.get("dt", dt);
  parser.get("time", time);
  parser.get("writeDirectory", writeDirectory);
  parser.get("temperature", temperature);
  parser.get("boundary", boundary);

  if (!quiet && rank==0) {
    #if DEBUG==1
    cout << "Running in Debug mode.\n";
    #endif
    // Print SIMD type
    #if SIMD_TYPE==SIMD_NONE
    cout << "Not using SIMD.\n";
    #elif SIMD_TYPE==SIMD_SSE
    cout << "Using SSE.\n";
    #elif SIMD_TYPE==SIMD_AVX
    cout << "Using AVX.\n";
    #elif SIMD_TYPE==SIMD_AVX2
    cout << "Using AVX2.\n";
    #elif SIMD_TYPE==SIMD_MIC
    cout << "Using MIC.\n";
    #endif
  }
  
  // Record the time at which the program started.
  auto start_time = current_time();
  if (!quiet && rank==0) cout << "Starting up simulation...\n";

  // --- This creator creates gflow simulations
  Creator *creator = nullptr;
  // Assign a specific type of creator
  if (bipartite_flag)  creator = new BipartiteBoxCreator(&parser);
  else if (debug_flag) creator = new DebugCreator(&parser);
  else if (load!="") {
    creator = new FileParseCreator(&parser, load);
  }
  else creator = new BoxCreator(&parser);

  // Set dimensions
  creator->setDimensions(dimensions);

  // --- Set boundary conditions
  switch (boundary) {
    case 0: {
      creator->setBCFlag(BCFlag::OPEN);
      break;
    }
    default:
    case 1: {
      creator->setBCFlag(BCFlag::WRAP);
      break;
    }
    case 2: {
      creator->setBCFlag(BCFlag::REFL);
      break;
    }
    case 3: {
      creator->setBCFlag(BCFlag::REPL);
      break;
    }
  }

  // --- Create a gflow simulation
  GFlow *gflow = creator->createSimulation();
  if (gflow==nullptr) {
    if (!quiet) cout << "GFlow was null. Exiting.\n";
    return 1;
  }

  // --- Make sure we didn't enter any illegal tokens - do this after gflow creation since creator uses flags
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    if (!quiet) cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }

  // --- Add data objects
  gflow->setStartRecTime(startRecTime);
  if (snapshot) gflow->addDataObject(new EndingSnapshot(gflow));
  if (totalKE || ke) gflow->addDataObject(new KineticEnergyData(gflow, ke));
  if (keTypes)     gflow->addDataObject(new KineticEnergyTypesData(gflow, true));
  if (bdForces)    gflow->addDataObject(new BoundaryForceData(gflow));
  if (timestep)    gflow->addDataObject(new TimeStepData(gflow));
  if (averages)    gflow->addDataObject(new AverageData(gflow));
  if (aveV) gflow->addDataObject(new AveVelocityData(gflow));
  if (aveP) gflow->addDataObject(new AveragePositionData(gflow));
  if (minDistances) gflow->addDataObject(new MinInteractingDistance(gflow));
  if (percolation) gflow->addDataObject(new PercolationData(gflow, skin));
  if (psnapshot) gflow->addDataObject(new PercolationSnapshot(gflow, skin));
  if (memdist)  gflow->addDataObject(new MemoryDistance(gflow));
  if (pressure) gflow->addDataObject(new PressureData(gflow));
  if (numberdata) gflow->addDataObject(new NumberData(gflow));
  // Add this last, as it takes the most time.
  if (animate) {
    auto pd = new PositionData(gflow);
    gflow->addDataObject(pd);
    if (videoLength>0) pd->setFPS((20.*videoLength)/time);
  }
  if (fps>0)    gflow->setFPS(fps); // Do after data objects are loaded
  gflow->setDMCmd(argc, argv);

  // --- Add modifiers
  if (temperature>0) gflow->addModifier(new TemperatureModifier(gflow, temperature));
  // Timestep adjustment
  if (target_steps>0) gflow->getIntegrator()->setTargetSteps(target_steps);
  if (step_delay>0) gflow->getIntegrator()->setStepDelay(step_delay);
  if (gravity!=0) gflow->addModifier(new ConstantAcceleration(gflow, gravity));
  if (damping) gflow->addModifier(new LinearVelocityDamping(gflow));

  // Set time step and request time
  if (!gflow->hasIntegrator()) {
    if (!quiet) cout << "GFlow does not have an integrator. Exiting.";
    return 0;
  }
  gflow->setDT(dt);
  gflow->requestTime(time);

  // Run the simulation
  if (!quiet && rank==0) {
    cout << "Initialized, ready to run:\t" << time_span(current_time(), start_time) << "\n";
    cout << "Running with " << gflow->getNumParticles() << " particles.\n";
  }
  if (gflow) {
    try {
      gflow->run();
    }
    catch (...) {
      if (!quiet && rank==0) cout << "Exited with exception.\n";
      // Write accumulated data to files
      if (rank==0) gflow->writeData(writeDirectory);
      // Rethrow the exception
      throw;
    }
    // More detailed exception handling
    // @todo Exception handling.
  }
  else {
    if (!quiet && rank==0) cout << "GFlow pointer was null. Exiting.\n";
    return 0;
  }
  if (!quiet && rank==0) cout << "Run is over:\t\t\t" << time_span(current_time(), start_time) << "\n";
  if (!quiet && rank==0) cout << "Ratio was:  \t\t\t" << gflow->getDataMaster()->getRatio() << "\n";

  // Write accumulated data to files
  if (rank==0) gflow->writeData(writeDirectory);
  if (!quiet && rank==0) cout << "Data write is over:\t\t" << time_span(current_time(), start_time) << "\n";

  // Delete creator, gflow
  if (creator) delete creator;
  if (gflow)   delete gflow;
  
  // Finalize mpi
  #if USE_MPI == 1
  #if _CLANG_ == 1
  MPI_Finalize();
  #else
  MPI::Finalize();
  #endif
  if (!quiet && rank==0) cout << "Finalized MPI.\n";
  #endif
  
  return 0;
}
