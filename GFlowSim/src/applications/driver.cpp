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
  #if DEBUG==1
  cout << "Running in Debug mode.\n";
  #endif
  
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
  
  // Record the time at which the program started.
  auto start_time = current_time();
  cout << "Starting up simulation...\n";

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
  cout << "Initialized MPI.\n";
  #endif

  // --- Options

  // Type of simulation
  bool bond_flag = false;
  bool bipartite_flag = false;
  bool debug_flag = false;
  bool flow_flag = false;
  bool crystal_flag = false;
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
  bool minDistances = false;
  bool percolation = false;
  bool psnapshot = false;
  bool memdist = false;
  bool pressure = false;
  bool numberdata = false;
  RealType skin = 0.;

  // Other options
  RealType gravity = 0.;
  bool adjustDT = false;
  RealType startRecTime = 0;
  RealType fps = -1.;
  RealType dt = 0.001;
  RealType time = 10.;
  string writeDirectory = "RunData";
  int boundary = 1;
  string monitor = ""; // Monitor file

  // Modifiers
  RealType temperature = 0;

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("bondbox", bond_flag);
  parser.get("bipartite", bipartite_flag);
  parser.get("debug", debug_flag); 
  parser.get("flow", flow_flag);
  parser.get("crystal", crystal_flag);
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
  parser.get("minDistances", minDistances);
  parser.get("percolation", percolation);
  parser.get("psnapshot", psnapshot);
  parser.get("memdist", memdist);
  parser.get("pressure", pressure);
  parser.get("numberdata", numberdata);
  parser.get("skin", skin);
  parser.get("gravity", gravity);
  parser.get("adjustDT", adjustDT);
  // parser.get("lj", adjustDT); // Adjust DT if lj is true
  parser.get("startRec", startRecTime);
  parser.get("fps", fps);
  parser.get("dt", dt);
  parser.get("time", time);
  parser.get("writeDirectory", writeDirectory);
  parser.get("temperature", temperature);
  parser.get("boundary", boundary);

  // --- This creator creates gflow simulations
  Creator *creator = nullptr;
  // Assign a specific type of creator
  if      (bond_flag)      creator = new BondBoxCreator(&parser);
  else if (bipartite_flag) creator = new BipartiteBoxCreator(&parser);
  else if (debug_flag)     creator = new DebugCreator(&parser);
  else if (flow_flag)      creator = new FlowCreator(&parser);
  else if (crystal_flag)   creator = new LineCrystalCreator(&parser);
  else if (load!="") {
    creator = new FileParseCreator(&parser, load);
  }
  else                     creator = new BoxCreator(&parser);

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
    cout << "GFlow was null. Exiting.\n";
    return 1;
  }

  // --- Make sure we didn't enter any illegal tokens - do this after gflow creation since creator uses flags
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }

  // --- Add data objects
  gflow->setStartRecTime(startRecTime);
  if (snapshot) gflow->addDataObject(new EndingSnapshot(gflow));
  if (sectorData)  gflow->addDataObject(new SectorizationData(gflow));
  if (totalKE || ke) gflow->addDataObject(new KineticEnergyData(gflow, ke));
  if (keTypes)     gflow->addDataObject(new KineticEnergyTypesData(gflow, true));
  if (secRemake)   gflow->addDataObject(new SectorizationRemakeData(gflow));
  if (bdForces)    gflow->addDataObject(new BoundaryForceData(gflow));
  if (timestep)    gflow->addDataObject(new TimeStepData(gflow));
  if (averages)    gflow->addDataObject(new AverageData(gflow));
  if (minDistances) gflow->addDataObject(new MinInteractingDistance(gflow));
  if (percolation) gflow->addDataObject(new PercolationData(gflow, skin));
  if (psnapshot) gflow->addDataObject(new PercolationSnapshot(gflow, skin));
  if (memdist)  gflow->addDataObject(new MemoryDistance(gflow));
  if (pressure) gflow->addDataObject(new PressureData(gflow));
  if (numberdata) gflow->addDataObject(new NumberData(gflow));
  // Add this last, as it takes the most time.
  if (animate)  gflow->addDataObject(new PositionData(gflow));
  if (fps>0)    gflow->setFPS(fps); // Do after data objects are loaded
  gflow->setDMCmd(argc, argv);

  // --- Add modifiers
  if (temperature>0) gflow->addModifier(new TemperatureModifier(gflow, temperature));
  // Timestep adjustment
  if (adjustDT) gflow->addModifier(new TimestepModifier(gflow));
  if (gravity!=0) gflow->addModifier(new ConstantAcceleration(gflow, gravity));

  // Set time step and request time
  if (!gflow->hasIntegrator()) {
    cout << "GFlow does not have an integrator. Exiting.";
    return 0;
  }
  gflow->setDT(dt);
  gflow->requestTime(time);

  // Run the simulation
  cout << "Initialized, ready to run:\t" << time_span(current_time(), start_time) << "\n";
  if (gflow) gflow->run();
  else {
    cout << "GFlow pointer was null. Exiting.\n";
    return 0;
  }
  cout << "Run is over:\t\t\t" << time_span(current_time(), start_time) << "\n";
  cout << "Ratio was:  \t\t\t" << gflow->getDataMaster()->getRatio() << "\n";

  // Write accumulated data to files
  gflow->writeData(writeDirectory);
  cout << "Data write is over:\t\t" << time_span(current_time(), start_time) << "\n";

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
  cout << "Finalized MPI.\n";
  #endif

  return 0;
}
