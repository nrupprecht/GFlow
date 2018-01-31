#include "FractureSimulation.hpp"

namespace GFlow {

  FractureSimulation::FractureSimulation() : SimulationBase(), runTime(10), delay(0.01), markBreaks(true) {};

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
    parser.get("markBreaks", markBreaks);

    // Standard parsing options and checking
    standardParsing();
  }

  // Run the simulation
  void FractureSimulation::run() {
    if (integrator==nullptr) return;

    // Integration
    bool running = true;
    RealType time = 0, timer = 0;

    // Initialize
    integrator->initialize(runTime);

    // Set up verlet list tracking (do after integrator initialization)
    int num = simData->getDomainSize();
    bondLists.initialize(simData);
    RealType cut = bondLists.getCutoff();
    bondLists.setCutoff(0.5*cut); // Set the cutoff and skin depth
    bondLists.createVerletLists(); // Create the initial verlet lists
    verletListRecord.resize(num);
    accessor.resize(num);
    setVerletListRecord();

    // Pre
    integrator->preIntegrate();

    // Run the simulation
    while (running) {
      // Perform an integration step
      integrator->step();
      // Update the time and running flag
      time = integrator->getTime();
      running = integrator->isRunning();
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
    if (integrator==nullptr) return;
    if (integrator->getSectorization()==nullptr) return;
    // Update the verlet lists and get them
    bondLists.createVerletLists();
    auto& verletList = bondLists.getVerletList();

    // Clear the marked set
    marked.clear();

    // Set up the accesser
    for (auto &a : accessor) a = nullptr;
    for (auto &vl : verletList) accessor.at(*vl.begin()) = &vl;

    // Keep track of new bondings and breaking bonds
    vector<pair<int, int> > bonds, breaks;

    // Look for bonds and breaks
    compare(bonds, breaks);

    // Set verlet list record to be the current verlet list. We have to sort the data
    setVerletListRecord(verletList);

    // Get the current time, for recording when this happened
    RealType currentTime = integrator->getTime();
    // If there is bond data to record
    if (!bonds.empty()) {
      // Add bond pairs
      for (auto &b : bonds) {
	bondings.push_back(std::tuple<RealType,int,int>(currentTime, b.first, b.second));
	// We keep an unordered set of all the particles bonding or undergoing breakages
	if (!markBreaks) {
	  marked.insert(b.first);
	  marked.insert(b.second);
	}

      }
      // Update cumulative counting statistics
      numBonds += bonds.size();
      cumNumBonds.push_back(pair<RealType,int>(currentTime, numBonds));
    }
    // If there is breakage data to record
    if (!breaks.empty()) {
      // Add break pairs
      for (auto &b : breaks) {
	breakages.push_back(std::tuple<RealType,int,int>(currentTime, b.first, b.second));
	// We keep an unordered set of all the particles bonding or undergoing breakages
	if (markBreaks) {
	  marked.insert(b.first);
	  marked.insert(b.second);
	}
      }
      // Update cumulative counting statistics
      numBreaks += breaks.size();
      cumNumBreaks.push_back(pair<RealType,int>(currentTime, numBreaks));
    }

    // Give the data record the marked particles
    dataRecord->setMarked(marked);
  }

  // Checks for bonds and breaks by comparing current and past verlet lists
  inline void FractureSimulation::compare(vector<pair<int,int> >& bonds, vector<pair<int,int> >& breaks) {
    // Lambda to find an element in a list
    auto check = [] (const int n, const VListSubType *list)->bool {
      if (list==nullptr) return false;
      for (auto &m : *list) if (m==n) return true;
      return false;
    };

    // Look for breaks
    for (int i=0; i<verletListRecord.size(); ++i) {
      if (!verletListRecord.at(i).empty()) {
	// For all of [i]'s neighbors, check that they stay neighbors
	for (auto nb : verletListRecord.at(i)) {
	  // Doesn't matter if you are in your verlet list. This is not a bonding/breakage
	  if (nb==i) continue;
	  // Check if [nb] is still in [i]'s neighbor list this time, or if [i] is in [nb]'s list this time
	  if (!check(nb, accessor.at(i)) && 
	      !check(i, accessor.at(nb)))
	    breaks.push_back(pair<int,int>(i,nb));
	}
      }
    }

    // Look for bonds
    for (int i=0; i<accessor.size(); ++i) {
      if (accessor.at(i)!=nullptr) {
	// For all of [i]'s neighbors, check that they were already neighbors. If not, a bond has formed
	for (auto nb : *accessor.at(i))
	  // Check if [nb] was in [i]'s neighbor list last time, or if [i] was in [nb]'s neighbor list last time
	  if (!check(nb, &verletListRecord.at(i)) && 
	      !check(i, &verletListRecord.at(nb)))
	    bonds.push_back(pair<int,int>(i,nb));
      }
    }
  }

  inline void FractureSimulation::setVerletListRecord(VListType& verletList) {
    // Clear the old data
    for (int i=0; i<verletListRecord.size(); ++i) verletListRecord.at(i).clear();
    // Clear accessor
    for (int i=0; i<accessor.size(); ++i) accessor.at(i) = nullptr;
    // Insert new data
    for (const auto &nl : verletList) {
      int head = *nl.begin();
      verletListRecord.at(head) = nl; // If this causes an error, then there are new particles
    }
  }

  inline void FractureSimulation::setVerletListRecord() {
    // Do checks
    if (integrator==nullptr) return;
    if (integrator->getSectorization()==nullptr) return;
    // Set up verlet list recording
    auto& verletList = bondLists.getVerletList();
    // Set the verlet list record
    setVerletListRecord(verletList);
  }

}
