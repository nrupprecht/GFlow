// For argument parsing
#include "utility/ArgParse.hpp"

// Simulation creators
#include "allcreators.hpp"

// All data objects to choose from
#include "alldataobjects.hpp"

// All the modifiers
#include "allmodifiers.hpp"

// MPI communicator.
#include "parallel/mpi-communication.hpp"

// For strlen
#include <cstring>

// For std::bad_alloc
#include <new>

/*
*  --- NOTES:
*  Running functions dynamically:
*  <https://stackoverflow.com/questions/11016078/is-it-possible-to-create-a-function-dynamically-during-runtime-in-c>
*  <http://burnttoys.blogspot.com/2011/04/how-to-allocate-executable-memory-on.html>
*
*  Using XCode: <https://hiltmon.com/blog/2015/08/01/simple-c-plus-plus-from-makefiles-to-xcode-builds/>
*
*/

using namespace GFlowSimulation;

#include "other/evaluation.hpp"

void finalize_mpi() {
  // Finalize mpi
  #if USE_MPI == 1
  #if _CLANG_ == 1
  MPI_Finalize();
  #else
  MPI::Finalize();
  #endif
  #endif
}

int main(int argc, char **argv) {
  // MPI related.
  int rank(0), numProc(1);

  // Allocate processor arguments. Do this outside USE_MPI.
  char **processor_argv = argv;

  #if USE_MPI == 1
  MPI_Init(&argc, &argv);
  // Get rank and number of processors.
  rank = MPIObject::getRank();
  numProc = MPIObject::getNumProc();

  // Make sure everyone gets the arguments supplied to the root processor.
  MPIObject::mpi_broadcast(argc);
  vector<int> arg_sizes(argc, 0);
  processor_argv = new char*[argc];

  // Broadcast to everyone else.
  if (rank==0) {
    // Get the sizes of each argument, so we can send them.
    for (int i=0; i<argc; ++i) arg_sizes[i] = strlen(argv[i]);
    // Broadcast sizes
    MPIObject::mpi_broadcast(arg_sizes);
    // Broadcast strings.
    for (int i=0; i<argc; ++i) {
      processor_argv[i] = new char[arg_sizes[i]+1];
      processor_argv[i][arg_sizes[i]] = '\0';
      copyVec(argv[i], processor_argv[i], arg_sizes[i]);
      MPIObject::mpi_broadcast(processor_argv[i], arg_sizes[i]);
    }
  }
  // Receive arguments from rank 0.
  else {
    // Receive sizes broadcast.
    MPIObject::mpi_broadcast(arg_sizes);
    // Receive strings.
    for (int i=0; i<argc; ++i) {
      processor_argv[i] = new char[arg_sizes[i]+1];
      processor_argv[i][arg_sizes[i]] = '\0';
      MPIObject::mpi_broadcast(processor_argv[i], arg_sizes[i]);
    }
  }

  // Wait here.
  MPIObject::barrier();
  #endif
  
  // --- Options

  // Type of simulation
  bool debug_flag = false;
  string load = "";

  // Data to gather
  bool animate = false; // Record positions
  real cavity = -1.f;
  bool snapshot = false;
  bool sectorData = false; // Record sector information
  bool ke = false; // Record kinetic energy
  bool rotE = false; // Record rotational kinetic energy
  real kebin = 0;
  bool energy = false; // Record total energy
  bool bondenergy = false;
  bool keTypes = false;
  bool totalKE = false; // Record average kinetic energy (per particle)
  bool aveOm = false; // Average angular velocity
  bool secRemake = false; 
  bool bdForces = false;
  bool timestep = false;
  bool averages = false;
  bool aveV = false;
  bool aveP = false;
  bool dev = false;
  real cavity_stat = 0;
  bool minDistances = false;
  bool percolation = false;
  bool psnapshot = false;
  bool memdist = false;
  bool pressure = false;
  bool numberdata = false;
  bool phidata = false;
  bool stripex = false;
  bool centercorr = false;
  bool velocityvp = false;
  bool radiusvp = false;
  
  // Other options
  int dimensions = 2;
  real skin = 0.;
  bool quiet = false;
  real gravity = 0.;
  real skin_depth = -1.;
  bool damping = false;
  bool adjustDT = true;
  int target_steps = -1;
  int step_delay = -1;
  real startRecTime = 0;
  real fps = -1.;
  real videoLength = -1.;
  real dt = 0.001;
  real maxDT = -1;
  long double time = 10.;
  bool print = false;
  string writeDirectory = "RunData";
  int boundary = 1;
  bool forces = true;
  bool timing = true;
  string monitor = ""; // Monitor file
  real checkpoint = 0.f;
  // Modifiers
  real temperature = 0.;

  // <--- End of options.

  // Flag that we can use to sync whether any mpi rank encountered errors (not set by parser).
  bool no_errors = true;

  // --- For getting command line arguments
  ArgParse parser;
  try {
    parser.set(argc, processor_argv);
  }
  catch (ArgParse::IllegalToken tok) {
    cout << "Illegal token in command line argument: " << tok.str << ". Exiting.\n";
    no_errors = false;
  }
  // Make sure this gets deleted even if there was an error and processor_argv is dynamically allocated.
  #if USE_MPI == 1
    for (int i=0; i<argc; ++i) delete [] processor_argv[i];
    delete [] processor_argv;
  #endif // USE_MPI == 1
  // Now that that's done, check for errors.
  MPIObject::mpi_and(no_errors);
  if (!no_errors) {
    finalize_mpi();
    return 1;
  }
  
  parser.get("debug", debug_flag); 
  parser.get("load", load);
  parser.get("animate", animate);
  parser.get("cavity", cavity);
  parser.get("snapshot", snapshot);
  parser.get("sectorData", sectorData);
  parser.get("KE", ke);
  parser.get("RKE", rotE);
  parser.get("KEBin", kebin);
  parser.get("energy", energy);
  parser.get("bondenergy", bondenergy);
  parser.get("KETypes", keTypes);
  parser.get("totalKE", totalKE);
  parser.get("aveOm", aveOm);
  parser.get("secRemake", secRemake);
  parser.get("bdForces", bdForces);
  parser.get("timestep", timestep);
  parser.get("averages", averages);
  parser.get("aveV", aveV);
  parser.get("aveP", aveP);
  parser.get("dev", dev);
  parser.get("cavity-stat", cavity_stat);
  parser.get("minDistances", minDistances);
  parser.get("percolation", percolation);
  parser.get("psnapshot", psnapshot);
  parser.get("memdist", memdist);
  parser.get("pressure", pressure);
  parser.get("numberdata", numberdata);
  parser.get("phidata", phidata);
  parser.get("stripex", stripex);
  parser.get("centercorr", centercorr);
  parser.get("velocityvp", velocityvp);
  parser.get("radiusvp", radiusvp);
  parser.get("dimensions", dimensions);
  parser.get("skin", skin);
  parser.get("quiet", quiet);
  parser.get("gravity", gravity);
  parser.get("skindepth", skin_depth);
  parser.get("damping", damping);
  parser.get("adjustDT", adjustDT);
  parser.get("target_steps", target_steps);
  parser.get("step_delay", step_delay);
  parser.get("startRec", startRecTime);
  parser.get("fps", fps);
  parser.get("videoLength", videoLength);
  parser.get("dt", dt);
  parser.get("maxDT", maxDT);
  parser.get("time", time);
  parser.get("print", print);
  parser.get("writeDirectory", writeDirectory);
  parser.get("temperature", temperature);
  parser.get("boundary", boundary);
  parser.get("forces", forces);
  parser.get("timing", timing);
  parser.get("checkpoint", checkpoint);

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

  // Seed the global generator.
  seedGlobalGenerator();

  // --- This creator creates gflow simulations
  Creator *creator = nullptr;
  // Assign a specific type of creator
  if (debug_flag) creator = new DebugCreator(&parser);
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
  GFlow *gflow = nullptr;
  try {
    gflow = creator->createSimulation();
  }
  catch (const Exception exc) {
    cout << "Simulation Creation: Exception caught on rank: " << rank << ", Message: " << exc.message << "\n";
    no_errors = false;
  }
  catch (...) {
    cout << "Simulation Creation: Creator encountered error (not inheriting from exception) on rank " << rank << " while creating the simulation. Exiting.\n";
    no_errors = false;
  }
  if (gflow==nullptr && no_errors) {
    if (!quiet) cout << "Simulation Creation: GFlow was null on rank " << rank << ". Exiting.\n";
    no_errors = false;
  }
  MPIObject::mpi_and(no_errors);
  if (!no_errors) {
    finalize_mpi();
    return 1;
  }

  // --- Make sure we didn't enter any illegal tokens - do this after gflow creation since creator uses flags

  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    // Since an illegal option will be illegal on every processor, we only have rank 0 print out an error message.
    if (!quiet && rank==0) cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    no_errors = false;
  }
  MPIObject::mpi_and(no_errors);
  if (!no_errors) {
    finalize_mpi();
    return 1;
  }

  // Whether we want to do full timing or not.
  TimedObject::setTimingOn(timing);

  // --- Add data objects
  gflow->setStartRecTime(startRecTime);
  if (snapshot)    gflow->addDataObject(make_shared<EndingSnapshot>(gflow));
  if (totalKE||ke) gflow->addDataObject(make_shared<KineticEnergyData>(gflow, ke));
  if (rotE)        gflow->addDataObject(make_shared<RotationalEnergyData>(gflow));
  if (kebin>0)     gflow->addDataObject(make_shared<KineticEnergyBin>(gflow, kebin));
  if (energy)      gflow->addDataObject(make_shared<TotalEnergyData>(gflow));
  if (bondenergy)  gflow->addDataObject(make_shared<BondedEnergyData>(gflow));
  if (keTypes)     gflow->addDataObject(make_shared<KineticEnergyTypesData>(gflow, true));
  if (aveOm)       gflow->addDataObject(make_shared<AverageOmegaData>(gflow));
  if (bdForces)    gflow->addDataObject(make_shared<BoundaryForceData>(gflow));
  if (timestep)    gflow->addDataObject(make_shared<TimeStepData>(gflow));
  if (averages)    gflow->addDataObject(make_shared<AverageData>(gflow));
  if (aveV)        gflow->addDataObject(make_shared<AverageVelocityData>(gflow));
  if (aveP)        gflow->addDataObject(make_shared<AveragePositionData>(gflow));
  if (dev)         gflow->addDataObject(make_shared<OscillationData>(gflow));
  if (0<cavity_stat) gflow->addDataObject(make_shared<CavityStatistics>(gflow, cavity_stat));
  if (minDistances)gflow->addDataObject(make_shared<MinInteractingDistance>(gflow));
  if (percolation) gflow->addDataObject(make_shared<PercolationData>(gflow, skin));
  if (psnapshot)   gflow->addDataObject(make_shared<PercolationSnapshot>(gflow, skin));
  if (memdist)     gflow->addDataObject(make_shared<MemoryDistance>(gflow));
  if (pressure)    gflow->addDataObject(make_shared<PressureData>(gflow));
  if (numberdata)  gflow->addDataObject(make_shared<NumberData>(gflow));
  if (phidata)     gflow->addDataObject(make_shared<PhiData>(gflow));
  if (centercorr)  gflow->addDataObject(make_shared<CenterCorrelation>(gflow));
  if (velocityvp)  gflow->addDataObject(make_shared<VelocityVolumePlot>(gflow));
  if (radiusvp)    gflow->addDataObject(make_shared<RadiusVolumePlot>(gflow));
  if (stripex)     gflow->addModifier(make_shared<StripeX>(gflow));

  real videoTime = startRecTime>0 ? time - startRecTime : time;
  if (videoLength<0) videoLength = videoTime;
  if (animate) {
    auto pd = make_shared<PositionData>(gflow);
    gflow->addDataObject(pd);
    if (stripex) {
      pd->add_scalar_data_entry("StripeX");
    }
    if (videoTime>0) pd->setFPS((20.*videoLength)/videoTime);
  }
  if (0<cavity) {
    auto cv = make_shared<CavityPositions>(gflow, cavity);
    gflow->addDataObject(cv);
    if (videoTime>0) cv->setFPS((20.*videoLength)/videoTime);
  }
  if (fps>0) gflow->setFPS(fps); // Do after data objects are loaded
  gflow->setDMCmd(argc, argv);
  gflow->setPrintUpdates(print);

  // --- Add modifiers
  if (temperature>0) gflow->addModifier(make_shared<TemperatureModifier>(gflow, temperature));
  // Add or modify some things.
  if (gravity!=0)     gflow->addModifier(make_shared<ConstantAcceleration>(gflow, gravity));
  if (damping)        gflow->addModifier(make_shared<LinearVelocityDamping>(gflow));
  if (skin_depth>0)   gflow->getInteractionHandler()->setSkinDepth(skin_depth);
  if (step_delay>0)   gflow->getIntegrator()->setStepDelay(step_delay);
  if (target_steps>0) gflow->getIntegrator()->setTargetSteps(target_steps);

  // Turns off forces.
  if (!forces) gflow->setUseForces(false);

  // Set time step and request time
  if (!gflow->hasIntegrator()) {
    if (!quiet) cout << "GFlow does not have an integrator on rank " << rank << ". Exiting.";
    finalize_mpi();
    return 0;
  }
  gflow->setDT(dt);
  if (maxDT>0) gflow->setMaxDT(maxDT);
  gflow->getIntegrator()->setAdjustDT(adjustDT);
  gflow->requestTime(time);

  // Tell the data object what directory it is going to write to.
  if (checkpoint>0) {
    gflow->getDataMaster()->setDoCheckpoints(true);
    gflow->getDataMaster()->setCheckpointingDelay(checkpoint);
  }
  gflow->getDataMaster()->setWriteDirectory(writeDirectory);

  // Run the simulation
  int num_particles = gflow->getNumParticles();
  #if USE_MPI == 1
    // Sum all the particles on all the processors.
    MPIObject::mpi_sum0(num_particles);
  #endif
  if (!quiet && rank==0) {
    cout << "Initialized, ready to run:\t" << time_span(current_time(), start_time) << "\n";
    cout << "Running with " << num_particles << " particles.\n";
  }
  try {
    gflow->run();
  }
  catch (const std::bad_alloc& err) {
    cout << "Exception caught on rank: " << rank << ": std::bad_alloc. Message: " << err.what() << "\n";
    no_errors = false;
  }
  catch (const Exception exc) {
    cout << "Exception caught on rank: " << rank << ", Message: " << exc.message << "\n";
    no_errors = false;
  }
  catch (...) {
    cout << "Exception (not inheriting from Exception class) caught on rank.\n";
    no_errors = false;
  }

  MPIObject::mpi_and(no_errors);
  if (!no_errors) {
    cout << "One or more processors exited with exception. Exiting.\n";
    // Write accumulated data to files
    if (rank==0) gflow->writeData(writeDirectory);
    finalize_mpi();
    return 1;
  }

  #if USE_MPI == 1
  // Sync here.
  MPIObject::barrier();
  #endif // USE_MPI == 1

  if (!quiet && rank==0) cout << "Run is over:\t\t\t" << time_span(current_time(), start_time) << "\n";
  if (!quiet && rank==0) cout << "Ratio was:  \t\t\t" << gflow->getDataMaster()->getRatio() << "\n";

  // Write accumulated data to files.
  // Only rank 0 will actually write anything, but information will need to sync, so all processes must call this function.
  gflow->writeData(writeDirectory); 

  // Print message that data write is over.
  if (!quiet && rank==0) cout << "Data write is over:\t\t" << time_span(current_time(), start_time) << "\n";

  // Delete creator, gflow
  if (creator) delete creator;
  if (gflow)   delete gflow;

  // Final barrier, to make sure that all data has been written.
  #if USE_MPI == 1
  // Sync here.
  MPIObject::barrier();
  finalize_mpi();
  #endif // USE_MPI == 1
  
  return 0;
}
