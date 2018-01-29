#include "FractureSimulation.hpp"

namespace GFlow {

  FractureSimulation::FractureSimulation() : SimulationBase(), runTime(10) {};

  void FractureSimulation::setUp(int argc, char** argv) {
    // Set the command line input
    SimulationBase::setUp(argc, argv);
    // Create data record, SimData, integrator
    dataRecord = new DataRecord;
    dataRecord->setCommand(argc, argv);
  }

  // Set parameters from the command line
  void FractureSimulation::parse() {
    // How long should the simulation run
    parser.get("time", runTime);

    // Standard parsing options and checking
    standardParsing();
  }

  // Run the simulation
  void FractureSimulation::run() {
    if (integrator==nullptr) return;

    // Initialize
    integrator->initialize(runTime);

    // Pre
    integrator->preIntegrate();

    // Integration
    bool running;
    RealType time = 0, timer = 0, delay = 0.1;

    // List of all neighbor lists for every particle. Two, one for current, one for new
    list<list<int> > neighborLists[2];
    int ol = 0, nl = 1;

    // Neighbor lists, indexed by head particle
    VListType neighborList; // ( vector<vector<int>> )

    while (running) {
      integrator->step();
      running = integrator->isRunning();
      time = integrator->getTime();

      // Check for fragmentation - just not every timestep
      if (time - timer > delay) {
	checkForBreaks();
	// Reset timer
	timer = time;
      }
    }

    // Post
    integrator->postIntegrate();
  }

  // Save data from the simulation
  void FractureSimulation::write() {
    standardWriting();
  }

  inline void FractureSimulation::checkForBreaks() {
    // Get a reference to the verlet lists
    auto& verletList = integrator->getSectorization()->getVerletList();

    // Keep track of new bondings and breaking bonds
    vector<pair<int, int> > bonds, breaks;

    // Compair to recorded verlet lists
    for (const auto &nl : verletList) {
      auto p = nl.begin();
      int head = *p; ++p;
      // Compair to the last recorded neighbor list for particle [head]
      auto list = verletListRecord.at(head);
      for (; p!=nl.end(); ++p) compare(list, nl, bonds, breaks);
    }

    // Set the verlet list record to be the current verlet list. We have to sort the data
    // Clear the old data
    for (auto &nl : verletListRecord) nl.clear();
    // Insert new data
    for (const auto &nl : verletList) {
      int head = *nl.begin();
      verletListRecord.at(head) = nl;
    }
    // Get the current time, for recording when this happened
    RealType currentTime = integrator->getTime();
    // If there is bond data to record
    if (!bonds.empty()) {
      // Add bond pairs
      for (auto &b : bonds) 
	bondings.push_back(std::tuple<RealType,int,int>(currentTime, b.first, b.second));
      // Update cumulative counting statistics
      numBonds += bonds.size();
      cumNumBonds.push_back(pair<RealType,int>(currentTime, numBonds));
    }
    // If there is breakage data to record
    if (!breaks.empty()) {
      // Add break pairs
      for (auto &b : breaks)
	breakages.push_back(std::tuple<RealType,int,int>(currentTime, b.first, b.second));
      // Update cumulative counting statistics
      numBreaks += breaks.size();
      cumNumBreaks.push_back(pair<RealType,int>(currentTime, numBreaks));
    }
  }

  inline void FractureSimulation::compare(const VListSubType& oldList, const VListSubType& newList, vector<pair<int,int> >& bonds, vector<pair<int,int> >& breaks) {
    
  }

}
