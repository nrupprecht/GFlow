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
  bool bipartite_flag = false;
  int nIters = 10;
  RealType width = 4.;
  RealType minPhi = 0.;
  RealType maxPhi = 0.8;
  RealType time = 60.;
  RealType startRec = 5.;
  RealType repulsion = 1.;

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("bipartite", bipartite_flag);
  parser.get("iters", nIters);
  parser.get("width", width);
  parser.get("minPhi", minPhi);
  parser.get("maxPhi", maxPhi);
  parser.get("time", time);
  parser.get("startRec", startRec);

  // Keep in bounds - max phi can be > 1 for bipartite
  minPhi = minPhi<0 ? 0 : minPhi;
  if (!bipartite_flag) maxPhi = maxPhi>1 ? 1 : maxPhi;

  RealType phi, invArea = 1./(2.*DIMENSIONS*pow(width,DIMENSIONS-1));
  vector<RPair> data;

  Creator *creator = nullptr;
  for (int iter=0; iter<nIters; ++iter) {
    // Set up
    phi = (maxPhi-minPhi)*static_cast<RealType>(iter)/nIters + minPhi;
    if (phi>0) {
      // --- Create a gflow simulation
      if (bipartite_flag) {
        BipartiteBoxCreator *bc = new BipartiteBoxCreator(&parser);
        bc->setWidth(width);
        bc->setPhi(phi);
        creator = bc;
      }
      else {
        BoxCreator *bc = new BoxCreator(&parser);
        bc->setWidth(width);
        bc->setPhi(phi);
        creator = bc;
      }
      GFlow *gflow = creator->createSimulation();
      if (gflow==nullptr) {
	cout << "GFlow was null. Exiting.\n";
	return 1;
      }
      gflow->setAllBCs(BCFlag::REPL);
      // Run
      BoundaryForceData *bfData = new BoundaryForceData(gflow);
      gflow->addDataObject(bfData);
      gflow->setRepulsion(repulsion*DEFAULT_HARD_SPHERE_REPULSION);
      gflow->run(time);
      // Data
      data.push_back(RPair(phi, invArea*bfData->getAverage()));
      // GFlow will delete bfData
      delete gflow, creator;
    }
    else data.push_back(RPair(0,0));
  }

  cout << "pr={";
  for (int i=0; i<data.size(); ++i) {
    cout << "{" << data.at(i).first << "," << data.at(i).second << "}";
    if (i!=data.size()-1) cout << ",";
  }
  cout << "};\n";

  return 0;
}
