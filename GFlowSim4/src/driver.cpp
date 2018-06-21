#include "boxcreator.hpp"

using namespace GFlowSimulation;

int main(int argc, char **argv) {

  // This creator creates gflow simulations
  Creator *creator = nullptr;

  // Assign a specific type of creator
  creator = new BoxCreator(argc, argv);

  // Create a gflow simulation
  GFlow *gflow = creator->createSimulation();

  // Run the simulation
  if (gflow) gflow->run(1.);
  else {
    cout << "GFlow pointer was null. Exiting.\n";
    return 0;
  }

  // Write accumulated data to files
  gflow->writeData("RunData");

  return 0;
}
