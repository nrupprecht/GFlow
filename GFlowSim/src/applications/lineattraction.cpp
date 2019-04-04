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
  RealType time = 10;
  RealType minh = 0;
  RealType maxh = 1;
  bool quiet = false;
  int iters = 10;
  int trials = 5;
  string writeDirectory = "LineAttraction/";
  // Other values
  string load = "configurations/lines.txt";

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("time", time);
  parser.get("minh", minh);
  parser.get("maxh", maxh);
  parser.get("quiet", quiet);
  parser.get("iters", iters);
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
  vector<RPair> data;

  // Create the directory
  rmdir(writeDirectory.c_str());
  mkdir(writeDirectory.c_str(), 0777);

  // Do many data gathering runs.
  RealType dh = (maxh - minh)/(iters-1);
  RealType h = minh;
  for (int i=0; i<iters; ++i) {
    // A creator.
    FileParseCreator creator(&parser, load);
    // Set variables in creator
    creator.setVariable("h", toStr(h)); // Line spacing
    // Accumulator for getting the average
    float ave = 0;
    // Run some trials
    for (int t=0; t<trials; ++t) {
      // Delete old gflow
      GFlow *gflow = creator.createSimulation();
      if (gflow==nullptr) throw false;

      // Add an ending snapshot object
      if (t==trials-1) gflow->addDataObject(new EndingSnapshot(gflow));

      auto pd = new PositionData(gflow);
      pd->clear_all_data_entries();
      pd->add_vector_data_entry("X");
      pd->add_scalar_data_entry("Sg");
      pd->add_integer_data_entry("Type");

      gflow->addDataObject(pd);
      RealType videoLength = 10;
      pd->setFPS((20.*videoLength)/time);

      // Find data objects.
      DataMaster *master = gflow->getDataMaster();
      auto &dataObjects = master->getDataObjects();
      // Find the group force object
      DataObject *dob = nullptr;
      for (auto obj : dataObjects) {
        if (obj->getName()=="LineEntropicForce") dob = obj;
      }
      if (dob==nullptr) throw false;

      // Run the program
      gflow->run(time);
      
      // Gather data
      auto *gnf = dynamic_cast<LineEntropicForce*>(dob);
      if (gnf==nullptr) throw false;
      ave += gnf->ave(0);

      // Write other data
      if (t==trials-1) gflow->writeData(writeDirectory+"data"+toStr(i));
      // Clean up
      delete gflow;
    }
    // Record
    ave /= static_cast<float>(trials);
    data.push_back(RPair(h, ave));

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

