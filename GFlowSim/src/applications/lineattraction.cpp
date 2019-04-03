// For argument parsing
#include "../utility/ArgParse.hpp"

// Simulation creators
#include "../allcreators.hpp"

// All data objects to choose from
#include "../alldataobjects.hpp"

// All the modifiers
#include "../allmodifiers.hpp"

using namespace GFlowSimulation;

int main(int argc, char **argv) {
  int rank(0);
  #if USE_MPI == 1
  int numProc(1);
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
  RealType time = 10;
  RealType minh = 0;
  RealType maxh = 1;
  bool quiet = false;
  int iters = 10;
  // Other values
  string load = "configurations/polymer.txt";

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("time", time);
  parser.get("minh", minh);
  parser.get("maxh", maxh);
  parser.get("quiet", quiet);
  parser.get("iters", iters);

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

  // --- Make sure we didn't enter any illegal tokens - do this after gflow creation since creator uses flags
  /*
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    if (!quiet) cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }
  */ 

  // Record data.
  vector<RPair> data;

  // Do many data gathering runs.
  RealType dh = (maxh - minh)/(iters-1);
  RealType h = minh;
  for (int i=0; i<iters; ++i) {
    // Set variables in creator
    FileParseCreator creator(&parser, load);
    creator.setVariable("h", toStr(h));
    creator.setVariable("pair", "1");
    // Create the simulation
    GFlow *gflow = creator.createSimulation();

    // If successful, run.
    if (gflow) {
      // Find data objects.
      DataMaster *master = gflow->getDataMaster();
      auto &dataObjects = master->getDataObjects();
      // Find the group force object
      DataObject *dob = nullptr;
      for (auto obj : dataObjects) {
        if (obj->getName()=="GroupForce") dob = obj;
      }

      // Run the program
      gflow->run(time);
      
      // Gather data
      auto *gnf = dynamic_cast<GroupNetForce*>(dob);
      if (gnf) data.push_back(RPair(h, gnf->ave(0)));
      // Delete gflow
      delete gflow;
    }
    // Increment h
    h += dh;
  }

  // Timing message
  cout << "Runs are over. Total time:\t\t\t" << time_span(current_time(), start_time) << "\n";

  // Print out data.
  cout << "data={";
  for (int i=0; i<data.size(); ++i) {
    cout << "{" << data[i].first << "," << data[i].second << "}";
    if (i!=data.size()-1) cout << ",";
  }
  cout << "};";
  
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

