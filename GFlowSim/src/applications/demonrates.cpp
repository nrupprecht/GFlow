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

inline void setDataObjects(GFlow *gflow, GraphObject* &KL, GraphObject* &KR, GraphObject* &NL, GraphObject* &NR, Parameters* &params) {
  // Get vector of data objects
  auto &dataObjects = gflow->getDataMaster()->getDataObjects();
  // Reset to null
  KL = KR = NL = NR = nullptr;
  params = nullptr;
  // Find the data object
  for (auto obj : dataObjects) {
    if      (obj->getName()=="KineticL") KL = dynamic_cast<GraphObject*>(obj);
    else if (obj->getName()=="KineticR") KR = dynamic_cast<GraphObject*>(obj);
    else if (obj->getName()=="NumberL" ) NL = dynamic_cast<GraphObject*>(obj);
    else if (obj->getName()=="NumberR" ) NR = dynamic_cast<GraphObject*>(obj);
    else if (obj->getName()=="Parameters") params = dynamic_cast<Parameters*>(obj);
  }
  if (KL==nullptr || KR==nullptr || NL==nullptr || NR==nullptr || params==nullptr) throw false;
}

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
  RealType time = 250;
  RealType phi = 0.05;
  RealType width = 40;
  RealType length = 0.5;
  RealType tau = 0.1;

  bool quiet = false;
  bool snapshot = false;
  int trials = 5;
  int bins = 25;
  string writeDirectory = "demon";
  // Other values
  string load = "configurations/demon.txt";

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("time", time);
  parser.get("phi", phi);
  parser.get("width", width);
  parser.get("length", length);
  parser.get("tau", tau);
  parser.get("quiet", quiet);
  parser.get("snapshot", snapshot);
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
  vector<RPair> average(bins, RPair(0,0));
  // All the data from the runs
  vector<vector<RealType> > allData;

  // Do many data gathering runs.
  
  // A creator.
  FileParseCreator creator(&parser, load);
  // Set variables in creator.
  creator.setVariable("tau", toStr(tau), true);
  creator.setVariable("width", toStr(width), true);
  creator.setVariable("length", toStr(length), true);
  // Pointers for data objects
  GraphObject *KL = nullptr, *KR = nullptr, *NL = nullptr, *NR = nullptr;
  Parameters *params = nullptr;
  

  // Iterate through different tau values.
  //for ()

  // Average energy and number currents.
  RealType aveE = 0, aveN = 0, stdE = 0, stdN = 0;
  vector<RealType> pointsE, pointsN;
  // Run some trials
  for (int t=0; t<trials; ++t) {
    // Delete old gflow
    GFlow *gflow = creator.createSimulation();
    if (gflow==nullptr) throw false;

    // Add an ending snapshot object
    if (snapshot && t==trials-1) gflow->addDataObject(new EndingSnapshot(gflow));

    // Find data objects.
    setDataObjects(gflow, KL, KR, NL, NR, params);

    // Run the program
    gflow->run(time);
    
    // Accumulate
    RealType aveDEDT = -1, aveDNDT = -1;
    if (KL->size()>0 && NL->size()>0) {
      aveDEDT = -(KL->last().second - KL->first().second)/(KL->last().first - KL->first().first);
      aveDNDT = -(NL->last().second - NL->first().second)/(NL->last().first - NL->first().first);
    }
    else throw false;
    // Update average accumulators and vectors.
    aveE += aveDEDT;
    aveN += aveDNDT;
    pointsE.push_back(aveE);
    pointsN.push_back(aveN);

    // Clean up
    delete gflow;
  }

  // Normalize
  aveE /= trials;
  aveN /= trials;

  // Compute sample standard deviation.

  for (auto en : pointsE) stdE += sqr(en - aveE);
  for (auto nm : pointsN) stdN += sqr(nm - aveN);
  stdE /= (trials-1);
  stdN /= (trials-1);
  stdE = sqrt(stdE);
  stdN = sqrt(stdN);

  cout << "Final: " << aveE << ", " << aveN << endl;
  cout << "Sample std dev: " << stdE/sqrt(trials) << ", " << stdN/sqrt(trials) << endl;

  // Timing message
  cout << "Runs are over. Total time:\t\t\t" << time_span(current_time(), start_time) << "\n";

  // Print average and std dev data
  ofstream fout(writeDirectory+"/forces.csv");
  
  fout.close();
  
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

