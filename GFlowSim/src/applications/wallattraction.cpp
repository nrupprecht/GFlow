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
  RealType time = 2000;
  RealType startRec = 200;
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
  vector<RPair> average(bins, RPair(0,0));
  // All the data from the runs
  vector<vector<RealType> > allData;

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
    vector<RealType> single_data(bins, 0);
    for (int i=0; i<bins; ++i) {
      average[i].first = entry[i].first;
      average[i].second += entry[i].second;
      single_data[i] = entry[i].second;
    }
    // Store data
    allData.push_back(single_data);

    // At the last trial
    if (t==trials-1) {
      // Create the directory
      mkdir(writeDirectory.c_str(), 0777);
      // Write data for the last trial
      gflow->writeData(writeDirectory+"/data"+toStr(t));
    }

    // Clean up
    delete gflow;
  }

  // Timing message
  cout << "Runs are over. Total time:\t\t\t" << time_span(current_time(), start_time) << "\n";

  // Make data into an average
  for (auto & datum : average) datum.second /= static_cast<double>(trials);
  // Find standard deviations
  vector<RealType> std(bins, 0);
  for (auto d : allData) {
    for (int i=0; i<bins; ++i) std[i] += sqr(d[i] - average[i].second);
  }
  for (int i=0; i<bins; ++i) std[i] = sqrt(std[i] / (trials - 1));

  // Print average and std dev data
  ofstream fout(writeDirectory+"/forces.csv");
  if (fout.fail()) {
    cout << "Ofstream failed to open file. Exiting.\n";
    return 0;
  }
  for (int i=0; i<average.size(); ++i) {
    fout << average[i].first << "," << average[i].second << "," << std[i] << endl;
  }
  // Close file stream.
  fout.close();

  // Print all data
  fout.open(writeDirectory+"/alldata.csv");
  if (fout.fail()) {
    cout << "Ofstream failed to open file. Exiting.\n";
    return 0;
  } 
  // First, print bins
  for (int i=0; i<bins; ++i) {
    fout << average[i].first;
    if (i!=bins-1) fout << ",";
  }
  fout << endl;
  // Then, print data (just <F>)
  for (const auto & datum : allData) {
    for (int i=0; i<datum.size(); ++i) {
      fout << datum[i];
      if (i!=datum.size()-1) fout << ",";
    }
    fout << endl;
  }
  // Close file stream.
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

