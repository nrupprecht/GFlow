// For argument parsing
#include "ArgParse.hpp"

// Simulation creators
#include "allcreators.hpp"

// All data objects to choose from
#include "alldataobjects.hpp"

// All the modifiers
#include "allmodifiers.hpp"

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
  bool binary_flag = false;
  bool debug_flag = false;

  // Data to gather
  bool animate; // Record positions
  bool vlData;  // Record verlet list information
  bool vlNums;  // Record numeric data about verlet lists
  bool sectorData; // Record sector information
  bool ke; // Record kinetic energy
  bool keTypes;
  bool aveKE; // Record average kinetic energy (per particle)
  bool secRemake; 
  RealType startRecTime = 0;
  RealType fps = 15.;
  RealType dt = 0.001;
  string writeDirectory = "RunData";

  // Modifiers
  RealType temperature = 0;

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("bondbox", bond_flag);
  parser.get("binarybox", binary_flag);
  parser.get("debug", debug_flag); 
  parser.get("animate", animate);
  parser.get("vlData", vlData);
  parser.get("vlNums", vlNums);
  parser.get("sectorData", sectorData);
  parser.get("KE", ke);
  parser.get("KETypes", keTypes);
  parser.get("aveKE", aveKE);
  parser.get("secRemake", secRemake);
  parser.get("startRec", startRecTime);
  parser.get("fps", fps);
  parser.get("dt", dt);
  parser.get("writeDirectory", writeDirectory);
  parser.get("temperature", temperature);

  // --- This creator creates gflow simulations
  Creator *creator = nullptr;
  // Assign a specific type of creator
  if (bond_flag)        creator = new BondBoxCreator(&parser);
  else if (binary_flag) creator = new BinaryBoxCreator(&parser);
  else if (debug_flag)  creator = new DebugCreator(&parser);
  else                  creator = new BoxCreator(&parser);

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
  if (animate) gflow->addDataObject(new PositionData(gflow));
  if (vlData)  gflow->addDataObject(new VerletListData(gflow));
  if (vlNums)  gflow->addDataObject(new VerletListNumberData(gflow));
  if (sectorData)  gflow->addDataObject(new SectorizationData(gflow));
  if (aveKE || ke) gflow->addDataObject(new KineticEnergyData(gflow, aveKE));
  if (keTypes)     gflow->addDataObject(new KineticEnergyTypesData(gflow, aveKE));
  if (secRemake)   gflow->addDataObject(new SectorizationRemakeData(gflow));
  gflow->setFPS(fps); // Do after data objects are loaded

  // --- Add modifiers
  if (temperature>0) gflow->addModifier(new TemperatureModifier(gflow, temperature));

  // Set time step
  gflow->setDT(dt);

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

  return 0;
}
