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
  RealType dt = 0.001;
  RealType time = 5.;
  RealType skin = 0.;
  bool adjustDT = false;

  // --- For getting command line arguments
  ArgParse parser(argc, argv);
  parser.get("bipartite", bipartite_flag);
  parser.get("iters", nIters);
  parser.get("width", width);
  parser.get("minPhi", minPhi);
  parser.get("maxPhi", maxPhi);
  parser.get("dt", dt);
  parser.get("time", time);
  parser.get("adjustDT", adjustDT);
  parser.get("lj", adjustDT); // If using lj, adjust dt

  // Keep in bounds - max phi can be > 1 for bipartite
  minPhi = minPhi<0 ? 0 : minPhi;
  if (!bipartite_flag) maxPhi = maxPhi>1 ? 1 : maxPhi;

  RealType phi, invArea = 1./(2.*DIMENSIONS*pow(width,DIMENSIONS-1));
  vector<RPair> data;

  Creator *creator = nullptr;
  for (int iter=0; iter<=nIters; ++iter) {
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
      // Repulsion force boundary conditions
      gflow->setAllBCs(BCFlag::REPL);
      // Run
      Clustering clustering(gflow);
      clustering.setSkin(skin);
      gflow->setDT(dt);
      // Timestep adjustment
      if (adjustDT) gflow->addModifier(new TimestepModifier(gflow));
      gflow->run(time);
      // Data
      clustering.findClusters();
      data.push_back(RPair(phi, static_cast<RealType>(clustering.getMaxClusterSize()) / gflow->getNumParticles()));
      // GFlow will delete bfData
      delete gflow, creator;
    }
    else data.push_back(RPair(0,0));
  }

  cout << "fraction={";
  for (int i=0; i<data.size(); ++i) {
    cout << "{" << data.at(i).first << "," << data.at(i).second << "}";
    if (i!=data.size()-1) cout << ",";
  }
  cout << "};\n";

  return 0;
}