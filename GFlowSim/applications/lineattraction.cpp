// For argument parsing
#include "utility/ArgParse.hpp"

// Simulation creators
#include "allcreators.hpp"

// All data objects to choose from
#include "alldataobjects.hpp"

// All the modifiers
#include "allmodifiers.hpp"

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
  RealType time = 60;
  RealType startRec = 0;
  RealType minh = 0;
  RealType maxh = 1.2;
  bool quiet = false;
  int iters = 13;
  int trials = 5;
  string writeDirectory = "LineAttraction";
  // Other values
  string load = "configurations/lines.txt";

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("time", time);
  parser.get("startRec", startRec);
  parser.get("minh", minh);
  parser.get("maxh", maxh);
  parser.get("quiet", quiet);
  parser.get("iters", iters);
  parser.get("trials", trials);
  parser.get("writeDirectory", writeDirectory);

  if (!quiet && rank == 0) {
    #if DEBUG == 1
    cout << "Running in Debug mode.\n";
    #endif
    // Print SIMD type
    #if SIMD_TYPE == SIMD_NONE
    cout << "Not using SIMD.\n";
    #elif SIMD_TYPE == SIMD_SSE
    cout << "Using SSE.\n";
    #elif SIMD_TYPE == SIMD_AVX
    cout << "Using AVX.\n";
    #elif SIMD_TYPE == SIMD_AVX2
    cout << "Using AVX2.\n";
    #elif SIMD_TYPE == SIMD_MIC
    cout << "Using MIC.\n";
    #endif
  }

  // Record the time at which the program started.
  auto start_time = current_time();

  // Record data.
  vector<vector<RealType> > data;

  // Create the directory
  //rmdir(writeDirectory.c_str());
  mkdir(writeDirectory.c_str(), 0777);

  // Do many data gathering runs.
  RealType dh = (maxh - minh) / (iters - 1);
  RealType h = minh;
  for (int i = 0; i < iters; ++i) {
    // A creator.
    FileParseCreator creator(&parser, load);
    // Set variables in creator
    creator.setVariable("h", toStr(h), true); // Line spacing
    // Accumulator for getting the average
    float ave = 0;
    vector<RealType> forces;
    // First entry is h
    forces.push_back(h);
    // Run some trials
    for (int t = 0; t < trials; ++t) {
      // Delete old gflow
      GFlow *gflow = creator.createSimulation();
      if (gflow == nullptr) {
        throw Exception("GFlow was null");
      }

      // Add an ending snapshot object
      if (t == trials - 1) {
        gflow->addDataObject(std::make_shared<EndingSnapshot>(gflow));
      }

      // Find data objects.
      DataMaster *master = gflow->getDataMaster();
      auto &dataObjects = master->getDataObjects();
      // Find the group force object
      std::shared_ptr<DataObject> dob;
      for (const auto& obj : dataObjects) {
        if (obj->getName() == "LineEntropicForce") {
          dob = obj;
        }
      }
      if (dob == nullptr) {
        throw Exception("Data object is null");
      }
      // Set start rec time.
      master->setStartRecTime(startRec);

      // Run the program
      gflow->run(time);

      // Gather data
      auto gnf = std::dynamic_pointer_cast<LineEntropicForce>(dob);
      if (gnf == nullptr) {
        throw Exception("Expected a Line entropic force object.");
      }
      RealType f = gnf->ave(0);;
      forces.push_back(f);

      // Write other data
      if (t == trials - 1) {
        gflow->writeData(writeDirectory + "/data" + toStr(i));
      }
      // Clean up
      delete gflow;
    }
    // Record
    data.push_back(forces);

    // Increment h
    h += dh;
  }

  // Timing message
  cout << "Runs are over. Total time:\t\t\t" << time_span(current_time(), start_time) << "\n";

  ofstream fout(writeDirectory + "/forces.csv");
  for (auto & i : data) {
    for (int j = 0; j < i.size(); ++j) {
      fout << i[j];
      if (j != i.size() - 1) {
        fout << ",";
      }
    }
    fout << endl;
  }

  // Finalize mpi
  #if USE_MPI == 1
  #if _CLANG_ == 1
  MPI_Finalize();
  #else
  MPI::Finalize();
  #endif
  if (!quiet && rank == 0) {
    cout << "Finalized MPI.\n";
  }
  #endif

  return 0;
}

