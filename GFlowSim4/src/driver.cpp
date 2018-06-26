// Simulation creators
#include "boxcreator.hpp"
#include "debugcreator.hpp"

// For argument parsing
#include "ArgParse.hpp"

// All data objects to choose from
#include "alldataobjects.hpp"

/*
*  --- NOTES:
*  Running functions dynamically:
*  <https://stackoverflow.com/questions/11016078/is-it-possible-to-create-a-function-dynamically-during-runtime-in-c>
*  <http://burnttoys.blogspot.com/2011/04/how-to-allocate-executable-memory-on.html>
*
*/

using namespace GFlowSimulation;

int main(int argc, char **argv) {

  // Type of simulation
  bool debug_flag = false;

  // Data to gather
  bool animate; // Record positions
  bool vlData;  // Record verlet list information

  // For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("debug", debug_flag);
  parser.get("animate", animate);
  parser.get("vlData", vlData);

  // This creator creates gflow simulations
  Creator *creator = nullptr;
  // Assign a specific type of creator
  if (debug_flag) creator = new DebugCreator(&parser);
  else            creator = new BoxCreator(&parser);

  // Create a gflow simulation
  GFlow *gflow = creator->createSimulation();

  // Make sure we didn't enter any illegal tokens - do this after gflow creation since creator uses flags
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }

  // Add data objects
  if (animate) gflow->addDataObject(new PositionData(gflow));
  if (vlData)  gflow->addDataObject(new VerletListData(gflow));

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
