// For argument parsing
#include "../utility/ArgParse.hpp"

// Simulation creators
#include "../allcreators.hpp"

// All data objects to choose from
#include "../alldataobjects.hpp"

// All the modifiers
#include "../allmodifiers.hpp"

#include <unistd.h>

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
  RealType time = 1000;
  RealType startRec = 100;
  RealType phi = 0.5;
  bool quiet = false;
  int trials = 5;
  int bins = 25;
  string writeDirectory = "WallAttraction";
  // Other values
  string load = "configurations/lines.txt";

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("time", time);
  parser.get("startRec", startRec);
  parser.get("phi", phi);
  parser.get("quiet", quiet);
  parser.get("trials", trials);
  parser.get("writeDirectory", writeDirectory);

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

  // Record data.
  vector<RPair> data(bins, RPair(0,0));

  // Create the directory
  //rmdir(writeDirectory.c_str());
  mkdir(writeDirectory.c_str(), 0777);

  // Do many data gathering runs.
  
  // A creator.
  FileParseCreator creator(&parser, load);
  // Set variables in creator.
  creator.setVariable("phi", toStr(phi), true); // Line density
  creator.setVariable("H",   toStr(2),   true); // Starting space
  
  // Run some trials
  for (int t=0; t<trials; ++t) {
    // Delete old gflow
    GFlow *gflow = creator.createSimulation();
    if (gflow==nullptr) throw false;

    // Add an ending snapshot object
    if (t==trials-1) gflow->addDataObject(new EndingSnapshot(gflow));

    // Find data objects.
    DataMaster *master = gflow->getDataMaster();
    auto &dataObjects = master->getDataObjects();
    // Find the data object
    DataObject *dob = nullptr;
    for (auto obj : dataObjects) {
      if (obj->getName()=="TwoWallBinForce") dob = obj;
    }
    if (dob==nullptr) throw false;
    // Pointer
    auto *gnf = dynamic_cast<TwoWallBinForce*>(dob);
    if (gnf==nullptr) throw false;
    gnf->setBins(bins);

    // Set start rec time.
    master->setStartRecTime(startRec);

    // Run the program
    gflow->run(time);
    
    // Accumulate
    auto entry = gnf->getEntry(0);
    cout << "{";
    for (int i=0; i<bins; ++i) {
      data[i].first = entry[i].first;
      data[i].second += entry[i].second;

      cout << "{" << data[i].first << "," << data[i].second << "}";
      if (i!=bins-1) cout << ",";
    }
    cout << "};\n";

    // Clean up
    delete gflow;
  }

  // Timing message
  cout << "Runs are over. Total time:\t\t\t" << time_span(current_time(), start_time) << "\n";

  ofstream fout(writeDirectory+"/forces.csv");
  for (int i=0; i<data.size(); ++i) {
    fout << data[i].first << "," << data[i].second / static_cast<double>(trials) << endl;
  }
  
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

