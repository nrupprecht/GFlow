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

  void FractureSimulation::checkForBreaks() {
    // Get a reference to the verlet lists
    auto& verletLists = integrator->getSectorization()->getVerletList();

    
  }

}
