// For argument parsing
#include "utility/ArgParse.hpp"
#include "allcreators.hpp"
#include "alldataobjects.hpp"
#include "allmodifiers.hpp"

#include <unistd.h>

// For tau, E, N, std E, std N
typedef std::tuple<float, float, float, float, float> datapoint;

typedef std::tuple<float, float, float, float, float, float, float> tuple7;

using namespace GFlowSimulation;

inline void setDataObjects(GFlow *gflow,
                           std::shared_ptr<GraphObject>& KL,
                           std::shared_ptr<GraphObject>& KR,
                           std::shared_ptr<GraphObject>& NL,
                           std::shared_ptr<GraphObject>& NR,
                           std::shared_ptr<GraphObject>& CE,
                           std::shared_ptr<GraphObject>& CN,
                           std::shared_ptr<Parameters>& params) {
  // Get vector of data objects
  const auto& dataObjects = gflow->getDataMaster()->getDataObjects();
  // Reset to null
  KL = KR = NL = NR = nullptr;
  params = nullptr;
  // Find the data object
  for (const auto &obj : dataObjects) {
    if (obj->getName() == "KineticL") {
      KL = std::dynamic_pointer_cast<GraphObject>(obj);
    }
    else if (obj->getName() == "KineticR") {
      KR = std::dynamic_pointer_cast<GraphObject>(obj);
    }
    else if (obj->getName() == "NumberL") {
      NL = std::dynamic_pointer_cast<GraphObject>(obj);
    }
    else if (obj->getName() == "NumberR") {
      NR = std::dynamic_pointer_cast<GraphObject>(obj);
    }
    else if (obj->getName() == "CurrentE") {
      CE = std::dynamic_pointer_cast<GraphObject>(obj);
    }
    else if (obj->getName() == "CurrentN") {
      CN = std::dynamic_pointer_cast<GraphObject>(obj);
    }
    else if (obj->getName() == "Parameters") {
      params = std::dynamic_pointer_cast<Parameters>(obj);
    }
  }
  if (KL == nullptr || KR == nullptr || NL == nullptr || NR == nullptr || params == nullptr) {
    throw Exception("One or more of the demon objects was null.");
  }
}

inline float kappa(float rho, float beta, float area, float mass, float tau) {
  return rho * tau * area / sqrt(2 * PI * beta * mass);
}

inline bool write_to_file(const vector<vector<pair<float, float> > > &all_data_type, const string &name) {
  std::ofstream fout(name);
  // If failure
  if (fout.fail()) {
    cout << "Failed to open file [" << name << "].\n";
    return false;
  }
  for (const auto &v : all_data_type) {
    // Write first
    for (const auto & i : v) {
      fout << i.first << ",";
    }
    // Write second
    for (int i = 0; i < v.size(); ++i) {
      fout << v[i].second;
      if (i != v.size() - 1) {
        fout << ",";
      }
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
  RealType time = 500;
  RealType width = 20;
  RealType length = 0.5;
  RealType tau = 0.1;
  int type = 0; // Demontype

  RealType minTau = 0.05;
  RealType maxTau = 2.0;
  int bins = 20;
  bool quiet = true;
  bool snapshot = false;
  int trials = 1;
  string writeDirectory = "demon";
  // Other values
  string load = "configurations/demon.txt";

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("time", time);
  parser.get("width", width);
  parser.get("length", length);
  parser.get("tau", tau);
  parser.get("type", type);
  parser.get("minTau", minTau);
  parser.get("maxTau", maxTau);
  parser.get("bins", bins);
  parser.get("quiet", quiet);
  parser.get("snapshot", snapshot);
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
  vector<datapoint> data;
  vector<tuple7> parameters;
  vector<pair<float, float> > prediction; // Predicted energy and number currents.
  vector<vector<pair<float, float> > > alldata; // Store each energy and number current for all runs.
  vector<vector<pair<float, float> > > allkappa; // Store all \kappa_r, \kappa_l
  vector<vector<pair<float, float> > > allratios; // Store all actual / prediction ratios.

  // Do many data gathering runs.

  // A creator.
  FileParseCreator creator(&parser, load);
  // Set variables in creator.
  creator.setVariable("width", toStr(width), true);
  creator.setVariable("length", toStr(length), true);
  creator.setVariable("type", toStr(type), true);

  // Pointers for data objects
  std::shared_ptr<GraphObject> KL, KR, NL, NR, CE, CN;
  std::shared_ptr<Parameters> params;

  RealType dt = 1 < bins ? (maxTau - minTau) / (bins - 1) : 0;
  // Iterate through different tau values.
  for (int i = 0; i < bins; ++i) {
    // Set tau
    RealType t = minTau + dt * i;
    creator.setVariable("tau", toStr(t), true);

    // Average energy and number currents.
    RealType aveE = 0, aveN = 0, stdE = 0, stdN = 0;
    vector<RealType> pointsE, pointsN;

    // Add a new entry to alldata
    alldata.emplace_back();
    allkappa.emplace_back();
    allratios.emplace_back();

    // Run some trials
    for (int tr = 0; tr < trials; ++tr) {
      // Delete old gflow
      auto gflow = creator.createSimulation();
      // Check to make sure nothing illegal was entered 
      parser.check();
      // Make sure gflow is not null
      if (gflow == nullptr) {
        throw Exception("GFlow object was null");
      }
      // Don't create any bmp images of plots. This makes it faster to download the data.
      // gflow->getDataMaster()->setAllPrintPlots(false);

      // Add an ending snapshot object
      if (snapshot && tr == trials - 1) {
        gflow->addDataObject(std::make_shared<EndingSnapshot>(gflow));
      }

      // Find data objects.
      setDataObjects(gflow, KL, KR, NL, NR, CE, CN, params);

      // Run the program
      gflow->run(time);

      // Accumulate
      if (KL->size() > 0 && NL->size() > 0) {
        RealType aveDEDT = CE->ave() / t;
        RealType aveDNDT = CN->ave() / t;
        // Update average accumulators and vectors.
        aveE += aveDEDT;
        aveN += aveDNDT;

        float area(0), mass(0), volume(0), energyL = KL->first().second, energyR = KR->first().second;
        float numL = NL->first().second, numR = NR->first().second;
        if (!params->find("Area", area)) {
          cout << "Couldn't find area.\n";
        }
        if (!params->find("Mass", mass)) {
          cout << "Couldn't find mass.\n";
        }
        if (!params->find("Vol", volume)) {
          cout << "Couldn't find volume.\n";
        }
        float rhoL = numL / volume, rhoR = numR / volume;
        float betaL = numL / energyL, betaR = numR / energyR;
        // Parameters
        parameters.emplace_back(rhoL, betaL, rhoR, betaR, t, area, mass);
        // Compute kappas
        float kappaL = kappa(rhoL, betaL, area, mass, t);
        float kappaR = kappa(rhoR, betaR, area, mass, t);
        allkappa.at(i).push_back(pair<float, float>(kappaL, kappaR));
        // Simple demon prediction
        float D = 3. / 2.;
        float I_energy = D * kappaL * exp(-kappaR) / (betaL * t);
        float I_number = I_energy * betaL / D;
        // Print a message
        if (!quiet) {
          cout << "Tau: " << t << " -- Ratios: " << aveDEDT / I_energy << ", " << aveDNDT / I_number << endl;
        }
        allratios.at(i).push_back(pair<float, float>(aveDEDT / I_energy, aveDNDT / I_number));
        pointsE.push_back(aveE);
        pointsN.push_back(aveN);
        alldata.at(i).push_back(pair<float, float>(aveDEDT, aveDNDT));
      }
      else {
        throw Exception();
      }

      if (i == bins - 1 && tr == trials - 1) {
        // Create the directory
        mkdir(writeDirectory.c_str(), 0777);
        // Write some data
        gflow->writeData(writeDirectory + "/data" + toStr(tr));
      }

      // Clean up
      delete gflow;
    }

    // Make sure directory has been created.
    mkdir(writeDirectory.c_str(), 0777);

    // Normalize
    if (!pointsE.empty()) {
      aveE /= pointsE.size();
    }
    if (!pointsN.empty()) {
      aveN /= pointsN.size();
    }

    // Compute sample standard deviation.
    for (auto en : pointsE) {
      stdE += sqr(en - aveE);
    }
    for (auto nm : pointsN) {
      stdN += sqr(nm - aveN);
    }
    stdE /= static_cast<RealType>(pointsE.size() - 1);
    stdN /= static_cast<RealType>(pointsN.size() - 1);
    stdE = sqrt(stdE);
    stdN = sqrt(stdN);

    // Print messages to screen
    if (!quiet) {
      cout << "Done with tau=" << t << endl;
      cout << "Final: " << aveE << ", " << aveN << endl;
      if (trials > 1) {
        cout << "Sample std dev: " << stdE / sqrt(trials) << ", " << stdN / sqrt(trials) << endl;
      }
      cout << endl;
    }

    // Store data
    data.emplace_back(t, aveE, aveN, stdE, stdN);
  }

  // Timing message
  cout << "Runs are over. Total time:\t\t\t" << time_span(current_time(), start_time) << "\n";

  // Print average and std dev data
  ofstream fout(writeDirectory + "/rates.csv");
  if (fout.fail()) {
    cout << "Failed to open file " << writeDirectory + "/rates.csv" << endl;
  }
  else {
    for (const auto& dp : data) {
      fout << std::get<0>(dp) << "," << std::get<1>(dp) << "," << std::get<2>(dp) << "," << std::get<3>(dp) << ","
           << std::get<4>(dp) << "\n";
    }
  }
  fout.close();

  // Print average and std dev data
  fout.open(writeDirectory + "/params.csv");
  if (fout.fail()) {
    cout << "Failed to open file " << writeDirectory + "/params.csv" << endl;
  }
  else {
    for (auto dp : parameters) {
      fout << std::get<0>(dp) << "," << std::get<1>(dp) << "," << std::get<2>(dp) << "," << std::get<3>(dp) << ","
           << std::get<4>(dp) << "," << std::get<5>(dp) << "," << std::get<6>(dp) << "\n";
    }
  }
  fout.close();

  write_to_file(alldata, writeDirectory + "/allrates.csv");
  write_to_file(allkappa, writeDirectory + "/allkappa.csv");
  write_to_file(allratios, writeDirectory + "/allratios.csv");

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

