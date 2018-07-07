// For argument parsing
#include "../src/ArgParse.hpp"

// Simulation creators
#include "../src/allcreators.hpp"

// All data objects to choose from
#include "../src/alldataobjects.hpp"

// All the modifiers
#include "../src/allmodifiers.hpp"

/*
*  --- NOTES:
*  Running functions dynamically:
*  <https://stackoverflow.com/questions/11016078/is-it-possible-to-create-a-function-dynamically-during-runtime-in-c>
*  <http://burnttoys.blogspot.com/2011/04/how-to-allocate-executable-memory-on.html>
*
*/

using namespace GFlowSimulation;

int main(int argc, char **argv) {
  // --- Options

  // Type of simulation
  bool bond_flag = false;
  bool bipartite_flag = false;
  bool debug_flag = false;

  // Data to gather
  bool animate = false; // Record positions
  bool snapshot = false;
  bool vlData = false;  // Record verlet list information
  bool vlNums = false;  // Record numeric data about verlet lists
  bool sectorData = false; // Record sector information
  bool ke = false; // Record kinetic energy
  bool keTypes = false;
  bool aveKE = false; // Record average kinetic energy (per particle)
  bool secRemake = false; 
  bool bdForces = false;
  RealType startRecTime = 0;
  RealType fps = 15.;
  RealType dt = 0.001;
  RealType time = 10.;
  string writeDirectory = "RunData";
  int boundary = 1;

  // Modifiers
  RealType temperature = 0;

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("bondbox", bond_flag);
  parser.get("bipartite", bipartite_flag);
  parser.get("debug", debug_flag); 
  parser.get("animate", animate);
  parser.get("snapshot", snapshot);
  parser.get("vlData", vlData);
  parser.get("vlNums", vlNums);
  parser.get("sectorData", sectorData);
  parser.get("KE", ke);
  parser.get("KETypes", keTypes);
  parser.get("aveKE", aveKE);
  parser.get("secRemake", secRemake);
  parser.get("bdForces", bdForces);
  parser.get("startRec", startRecTime);
  parser.get("fps", fps);
  parser.get("dt", dt);
  parser.get("time", time);
  parser.get("writeDirectory", writeDirectory);
  parser.get("temperature", temperature);
  parser.get("boundary", boundary);

  // --- This creator creates gflow simulations
  Creator *creator = nullptr;
  // Assign a specific type of creator
  if      (bond_flag)      creator = new BondBoxCreator(&parser);
  else if (bipartite_flag) creator = new BipartiteBoxCreator(&parser);
  else if (debug_flag)     creator = new DebugCreator(&parser);
  else                     creator = new BoxCreator(&parser);
  // --- Set boundary conditions
  switch (boundary) {
    case 0:
    creator->setBCFlag(BCFlag::OPEN);
    break;
    default:
    case 1:
    creator->setBCFlag(BCFlag::WRAP);
    break;
    case 2:
    creator->setBCFlag(BCFlag::REFL);
    break;
    case 3:
    creator->setBCFlag(BCFlag::REPL);
    break;
  }

  // --- Create a gflow simulation
  GFlow *gflow = creator->createSimulation();

  // --- Make sure we didn't enter any illegal tokens - do this after gflow creation since creator uses flags
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }

  // --- Add data objects
  gflow->setStartRecTime(startRecTime);
  if (animate)  gflow->addDataObject(new PositionData(gflow));
  if (snapshot) gflow->addDataObject(new EndingSnapshot(gflow));
  if (vlData)   gflow->addDataObject(new VerletListData(gflow));
  if (vlNums)   gflow->addDataObject(new VerletListNumberData(gflow));
  if (sectorData)  gflow->addDataObject(new SectorizationData(gflow));
  if (aveKE || ke) gflow->addDataObject(new KineticEnergyData(gflow, aveKE));
  if (keTypes)     gflow->addDataObject(new KineticEnergyTypesData(gflow, aveKE));
  if (secRemake)   gflow->addDataObject(new SectorizationRemakeData(gflow));
  if (bdForces)    gflow->addDataObject(new BoundaryForceData(gflow));
  gflow->setFPS(fps); // Do after data objects are loaded

  // --- Add modifiers
  if (temperature>0) gflow->addModifier(new TemperatureModifier(gflow, temperature));

  // Set time step and request time
  gflow->setDT(dt);
  gflow->requestTime(time);

  // Run the simulation
  if (gflow) gflow->run();
  else {
    cout << "GFlow pointer was null. Exiting.\n";
    return 0;
  }
  cout << "Run is over.\n";

  // Write accumulated data to files
  gflow->writeData(writeDirectory);
  cout << "Data write is over.\n";

  // Delete creator, gflow
  delete creator, gflow;

  return 0;
}
