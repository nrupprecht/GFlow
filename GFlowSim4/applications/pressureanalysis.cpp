// For argument parsing
#include "../src/ArgParse.hpp"

// Simulation creators
#include "../src/allcreators.hpp"

// All data objects to choose from
#include "../src/alldataobjects.hpp"

// All the modifiers
#include "../src/allmodifiers.hpp"

using namespace GFlowSimulation;

int main(int argc, char **argv) {
  // Options
  int nIters = 10;
  RealType width = 4.;
  RealType minPhi = 0.;
  RealType maxPhi = 0.8;
  RealType time = 60.;
  RealType startRec = 5.;

  // --- For getting command line arguments
  ArgParse parser(argc, argv);

  // --- Create a gflow simulation
  BoxCreator creator(&parser);
  creator.setWidth(width);

  RealType phi;
  vector<RPair> data;
  for (int iter=0; iter<nIters; ++iter) {
    // Set up
    phi = (maxPhi-minPhi)*static_cast<RealType>(iter)/nIters + minPhi;
    creator.setPhi(phi);
    GFlow *gflow = creator.createSimulation();
    // Run
    BoundaryForceData *bfData = new BoundaryForceData(gflow);
    gflow->addDataObject(bfData);
    gflow->run(time);
    // Data
    data.push_back(RPair(phi, bfData->getAverage()));
    // GFlow will delete bfData
    delete gflow;
  }

  return 0;
}
