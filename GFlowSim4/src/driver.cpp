#include "boxcreator.hpp"
#include "debugcreator.hpp"

/*
*
*  Running functions dynamically:
*  <https://stackoverflow.com/questions/11016078/is-it-possible-to-create-a-function-dynamically-during-runtime-in-c>
*  <http://burnttoys.blogspot.com/2011/04/how-to-allocate-executable-memory-on.html>
*
*/

using namespace GFlowSimulation;

int main(int argc, char **argv) {

  // This creator creates gflow simulations
  Creator *creator = nullptr;

  // Assign a specific type of creator
  // creator = new BoxCreator(argc, argv);
  creator = new DebugCreator(argc, argv);

  // Create a gflow simulation
  GFlow *gflow = creator->createSimulation();

  // Run the simulation
  if (gflow) gflow->run();
  else {
    cout << "GFlow pointer was null. Exiting.\n";
    return 0;
  }

  cout << "Run is over.\n";

  // Write accumulated data to files
  gflow->writeData("RunData");

  cout << "Data write is over.\n";

  return 0;
}
