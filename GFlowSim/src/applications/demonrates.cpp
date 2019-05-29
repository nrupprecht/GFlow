// For argument parsing
#include "../utility/ArgParse.hpp"

// Simulation creators
#include "../allcreators.hpp"

// All data objects to choose from
#include "../alldataobjects.hpp"

// All the modifiers
#include "../allmodifiers.hpp"

#include <unistd.h>

// For tau, E, N, std E, std N
typedef std::tuple<float, float, float, float, float> datapoint;

using namespace GFlowSimulation;

inline void setDataObjects(GFlow *gflow, GraphObject* &KL, GraphObject* &KR, GraphObject* &NL, GraphObject* &NR, GraphObject* &CE, GraphObject* &CN, Parameters* &params) {
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
    else if (obj->getName()=="CurrentE" ) CE = dynamic_cast<GraphObject*>(obj);
    else if (obj->getName()=="CurrentN" ) CN = dynamic_cast<GraphObject*>(obj);
    else if (obj->getName()=="Parameters") params = dynamic_cast<Parameters*>(obj);
  }
  if (KL==nullptr || KR==nullptr || NL==nullptr || NR==nullptr || params==nullptr) throw false;
}

inline float kappa(float rho, float beta, float area, float mass, float tau) {
  return rho * tau * area / sqrt(2 * PI * beta * mass);
}

inline bool write_to_file(const vector<vector<pair<float,float> > >& all_data_type, const string& name) {
  std::ofstream fout(name);
  // If failure
  if (fout.fail()) {
    cout << "Failed to open file [" << name << "].\n";
    return false;
  }
  for (const auto &v : all_data_type) {
    // Write first
    for (int i=0; i<v.size(); ++i) {
      fout << v[i].first << ",";
    }
    // Write second
    for (int i=0; i<v.size(); ++i) {
      fout << v[i].second;
      if (i!=v.size()-1) fout << ",";
    }
    fout << endl;
  }
  fout.close();
  // Return success.
  return true;
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

  RealType minTau = 0.05;
  RealType maxTau = 2.0;
  int bins = 25;
  bool quiet = true;
  bool snapshot = false;
  int trials = 5;
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
  parser.get("minTau", minTau);
  parser.get("maxTau", maxTau);
  parser.get("bins", bins);
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
  vector<datapoint> data;
  vector< vector<pair<float, float> > > alldata; // Store each energy and number current for all runs.
  vector< vector<pair<float, float> > > allkappa; // Store all \kappa_r, \kappa_l
  vector< vector<pair<float, float> > > allratios; // Store all actual / prediction ratios.

  // Do many data gathering runs.
  
  // A creator.
  FileParseCreator creator(&parser, load);
  // Set variables in creator.
  creator.setVariable("width", toStr(width), true);
  creator.setVariable("length", toStr(length), true);
  // Pointers for data objects
  GraphObject *KL = nullptr, *KR = nullptr, *NL = nullptr, *NR = nullptr, *CE = nullptr, *CN = nullptr;
  Parameters *params = nullptr;

  // Create the directory
  mkdir(writeDirectory.c_str(), 0777);
  
  RealType dt = bins>1 ? (maxTau - minTau)/(bins-1) : 0;
  // Iterate through different tau values.
  for (int i=0; i<bins; ++i) {
    // Set tau
    RealType t = minTau + dt*i;
    creator.setVariable("tau", toStr(t), true);

    // Average energy and number currents.
    RealType aveE = 0, aveN = 0, stdE = 0, stdN = 0;
    vector<RealType> pointsE, pointsN;

    // Add a new entry to alldata
    alldata.push_back( vector<pair<float, float> >());
    allkappa.push_back(vector<pair<float, float> >());
    allratios.push_back(vector<pair<float, float> >());

    // Run some trials
    for (int tr=0; tr<trials; ++tr) {
      // Delete old gflow
      GFlow *gflow = creator.createSimulation();
      if (gflow==nullptr) throw false;

      // Add an ending snapshot object
      if (snapshot && tr==trials-1) gflow->addDataObject(new EndingSnapshot(gflow));

      // Find data objects.
      setDataObjects(gflow, KL, KR, NL, NR, CE, CN, params);

      // Run the program
      gflow->run(time);
      
      // Accumulate
      if (KL->size()>0 && NL->size()>0) {
        RealType aveDEDT = CE->ave()/t; // -(KL->last().second - KL->first().second)/(KL->last().first - KL->first().first);
        RealType aveDNDT = CN->ave()/t; //-(NL->last().second - NL->first().second)/(NL->last().first - NL->first().first);
        // Update average accumulators and vectors.
        aveE += aveDEDT;
        aveN += aveDNDT;

        float area(0), mass(0), volume(0), energyL = KL->first().second, energyR = KR->first().second;
        float numL = NL->first().second, numR = NR->first().second;
        if (!params->find("Area", area)) cout << "Couldn't find area.\n";
        if (!params->find("Mass", mass)) cout << "Couldn't find mass.\n";
        if (!params->find("Vol", volume)) cout << "Couldn't find volume.\n";
        float rhoL = numL / volume, rhoR = numR / volume;
        float betaL = numL / energyL, betaR = numR / energyR;
        // Compute kappas
        float kappaL = kappa(rhoL, betaL, area, mass, t);
        float kappaR = kappa(rhoR, betaR, area, mass, t);
        allkappa.at(i).push_back(pair<float,float>(kappaL, kappaR));
        // Simple demon prediction
        float D = 3./2.;
        float I_energy = D * kappaL / t * exp(-kappaR) / betaL;
        float I_number = I_energy * betaL / D;

        cout << "Tau: " << t << " -- Ratios: " << aveDEDT / I_energy << ", " << aveDNDT / I_number << endl;

        allratios.at(i).push_back(pair<float,float>(aveDEDT / I_energy, aveDNDT / I_number));

        pointsE.push_back(aveE);
        pointsN.push_back(aveN);
        alldata.at(i).push_back(pair<float,float>(aveDEDT, aveDNDT));
      }
      else throw false;

      if (i==bins-1 && tr==trials-1) {
        // Write some data
        gflow->writeData(writeDirectory+"/data"+toStr(tr));
      }

      // Clean up
      delete gflow;
    }

    // Normalize
    if (!pointsE.empty()) aveE /= pointsE.size();
    if (!pointsN.empty()) aveN /= pointsN.size();

    // Compute sample standard deviation.
    for (auto en : pointsE) stdE += sqr(en - aveE);
    for (auto nm : pointsN) stdN += sqr(nm - aveN);
    stdE /= (pointsE.size()-1);
    stdN /= (pointsN.size()-1);
    stdE = sqrt(stdE);
    stdN = sqrt(stdN);

    // Print messages to screen
    if (!quiet) {
      cout << "Tau: " << t << endl;
      cout << "Final: " << aveE << ", " << aveN << endl;
      cout << "Sample std dev: " << stdE/sqrt(trials) << ", " << stdN/sqrt(trials) << endl;
      cout << endl;
    }

    // Store data
    data.push_back(datapoint(t, aveE, aveN, stdE, stdN));
  }

  // Timing message
  cout << "Runs are over. Total time:\t\t\t" << time_span(current_time(), start_time) << "\n";

  // Print average and std dev data
  ofstream fout(writeDirectory+"/rates.csv");
  if (fout.fail()) {
    cout << "Failed to open file " << writeDirectory+"/rates.csv" << endl;
  }
  else {
    for (auto dp : data)
      fout << std::get<0>(dp) << "," << std::get<1>(dp) << "," << std::get<2>(dp) << "," << std::get<3>(dp) << "," << std::get<4>(dp) << "\n";
    fout.close();
  }

  write_to_file(alldata, writeDirectory+"/allrates.csv");
  /*
  fout.open(writeDirectory+"/allrates.csv");
  if (fout.fail()) cout << "Failed to open file " << writeDirectory+"/allrates.csv" << endl;
  else {
    for (const auto &v : alldata) {
      // Write energy current
      for (int i=0; i<v.size(); ++i) {
        fout << v[i].first << ",";
      }
      // Write number current
      for (int i=0; i<v.size(); ++i) {
        fout << v[i].second;
        if (i!=v.size()-1) fout << ",";
      }
      fout << endl;
    }
  }
  */

  write_to_file(allkappa, writeDirectory+"/allkappa.csv");
  write_to_file(allratios, writeDirectory+"/allratios.csv");
  
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

