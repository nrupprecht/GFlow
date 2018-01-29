#include "StandardSimulation.hpp"

namespace GFlow {

  StandardSimulation::StandardSimulation() : SimulationBase(), runTime(10) {};
  
  void StandardSimulation::setUp(int argc, char** argv) {
    // Set the command line input
    SimulationBase::setUp(argc, argv);
    // Create data record, SimData, integrator
    dataRecord = new DataRecord;
    dataRecord->setCommand(argc, argv);
  }

  void StandardSimulation::parse() {
    // How long should the simulation run
    parser.get("time", runTime);
    // Standard parsing options and checking
    standardParsing();
  }

  void StandardSimulation::run() {
    integrator->integrate(runTime);
  }


  void StandardSimulation::write() {
    standardWriting();
  }
}
