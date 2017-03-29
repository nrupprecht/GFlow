/// driver.cpp --- Driver file
/// Nathaniel Rupprecht 2016
///

#include "../ArgParse.h"
#include "GFlowBase.h"
//#include "GFlow.h"

using MPI::COMM_WORLD;

#include <iostream>
using std::cout;
using std::endl;

int main(int argc, char** argv) {
  // Simulation parameters
  int number = -1;
  double width = 4;
  double height = 4;
  double radius = 0.05;
  double bR = 0.2;
  double density = 10;
  double velocity = 0.25;
  double dispersion = 0;
  double temperature = 0;
  double time = 1.;
  double start = 0;
  double phi = 0.5;
  double skinDepth = -1;
  int interaction = 0;
  bool interact = true;
  bool seedRand = true;

  // Animation Paramaters
  int mode     = 0;     // Animation mode
  bool animate = false;
  bool special = false;
  bool bubbles = false;
  bool visBubbles = false;
  bool omega   = false;
  bool KE      = false;
  bool LKE     = false;
  bool RKE     = false;
  bool cluster = false;
  bool triAlign = false;
  bool trackHeight  = false;
  bool novid   = false;
  bool printSectors = false;
  double fps = -1; // FPS

  // Simulation type
  bool square = true;
  bool buoyancy = false;

  // Load and save
  string loadFile = "";
  string loadBuoyancy = "";
  string saveFile = "";

  // Initialize MPI
  int rank, numProc;
  MPI::Init();
  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();
  MPI::COMM_WORLD.Barrier();

  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser;
  try {
    parser.set(argc, argv);
  }
  catch (ArgParse::IllegalToken token) {
    cout << "Illegal Token: " << token.c << ". Exiting.\n";
    exit(1);
  }
  parser.get("number", number);
  parser.get("width", width);
  parser.get("height", height);
  parser.get("radius", radius);
  parser.get("bR", bR);
  parser.get("density", density);
  parser.get("velocity", velocity);
  parser.get("dispersion", dispersion);
  parser.get("temperature", temperature);
  parser.get("time", time);
  parser.get("start", start);
  parser.get("phi", phi);
  parser.get("skinDepth", skinDepth);
  parser.get("interaction", interaction);
  parser.get("interact", interact);
  parser.get("srand", seedRand);
  parser.get("mode", mode);
  parser.get("animate", animate);
  parser.get("special", special);
  parser.get("bubbles", bubbles);
  parser.get("visBubbles", visBubbles);
  parser.get("omega", omega);
  parser.get("KE", KE);
  parser.get("LKE", LKE);
  parser.get("RKE", RKE);
  parser.get("cluster", cluster);
  parser.get("triAlign", triAlign);
  parser.get("trackHeight", trackHeight);
  parser.get("novid", novid);
  parser.get("printSectors", printSectors);
  parser.get("fps", fps);
  parser.get("square", square);
  parser.get("buoyancy", buoyancy);
  parser.get("loadFile", loadFile);
  parser.get("loadBuoyancy", loadBuoyancy);
  parser.get("saveFile", saveFile);
  // Make sure we didn't enter any illegal tokens (ones not listed above) on the command line
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }
  //----------------------------------------

  // Seed random number generators
  if (seedRand) {
    srand48( std::time(0) );
    srand( std::time(0) );
    seedNormalDistribution();
  }

  // Calculate number of particles given a packing fraction
  if (number==-1) {
    double Vol = width*height;
    number = Vol/(PI*sqr(radius))*phi;
  }

  /// Print condition summary
  if (rank==0) {
    cout << "----------------------- RUN SUMMARY -----------------------\n";
    cout << "Command: ";
    for (int i=0; i<argc; i++) cout << argv[i] << " ";
    cout << endl << endl; // Line break
    cout << "Using " << numProc << " processors.\n";
    cout << "  ..........................\n";
  }
  MPI::COMM_WORLD.Barrier();

  // Set up the simulation
  GFlowBase simulator;
  if (0<fps) simulator.setFrameRate(fps);
  if (0<skinDepth) simulator.setSkinDepth(skinDepth);
  // Create scenario
  if (loadFile=="" && loadBuoyancy=="") {
    if (buoyancy) simulator.createBuoyancyBox(radius, bR, density, width, height, velocity, dispersion, interaction);
    else if (square) simulator.createSquare(number, radius, width, height, velocity, dispersion, interaction);
    else throw false; // No selection
  }
  else if (loadBuoyancy!="") {
    if (!simulator.loadBuoyancy(loadBuoyancy, bR, velocity, density)) {
      cout << "Failed to load [" << loadFile << "], exiting.\n";
      return 0;
    }
  }
  else { // We must have that loadFile!=""
    if (!simulator.loadConfigurationFromFile(loadFile)) {
      cout << "Failed to load [" << loadFile <<"], exiting.\n";
      return 0;
    }
  }

  simulator.setTemperature(temperature);
  simulator.setDoInteractions(interact);
  simulator.setStartRec(start);

  simulator.setRecPositions(animate);
  simulator.setRecSpecial(special);
  simulator.setRecBubbles(bubbles);
  simulator.setVisBubbles(visBubbles);
  if (buoyancy || loadBuoyancy!="") simulator.setRestrictBubbleDomain(true);

  if (omega) simulator.addStatFunction(Stat_Omega, "omega");
  if (KE) simulator.addStatFunction(Stat_KE, "ke");
  if (LKE) simulator.addStatFunction(Stat_L_KE, "lke");
  if (RKE) simulator.addStatFunction(Stat_R_KE, "rke");
  if (cluster) simulator.addStatFunction(Stat_Clustering, "cluster");
  if (triAlign) simulator.addStatFunction(Stat_Triangle_Align, "triAlign");
  if (trackHeight) simulator.addStatFunction(Stat_Large_Object_Height, "height");

  // Get the actual number of particles in the simulation
  number = simulator.getSize();
  double filled = simulator.getFilledVolume(), vol = simulator.getVolume();
  if (rank==0) {
    cout << "Number: " << number << ", Packing fraction: " << filled/vol << endl;
    cout << "Dimensions: " << simulator.getWidth() << " x " << simulator.getHeight() << endl;
    cout << "Set up time: " << simulator.getSetUpTime() << endl;
    cout << "  ..........................\n";
  }

  if (interaction!=0) simulator.setInteractionType(interaction);

  // Run the simulation
  simulator.run(time);

  // Head node prints the run summary
  if (rank==0) {
    if (saveFile!="") {
      if (simulator.createConfigurationFile(saveFile)) cout << "Saved configuration to file [" << saveFile << "]" << endl;
    }

    double runTime = simulator.getRunTime(), transferTime = simulator.getTransferTime();
    int iters = simulator.getIter(), ndx = simulator.getNDX(), ndy = simulator.getNDY();
    int nsx = simulator.getNSX(), nsy = simulator.getNSY();
    cout << "Domains: " << ndx << " x " << ndy << ", Total: " << ndx*ndy << endl;
    cout << "Sectors: " << nsx << " x " << nsy << ", Per Domain: " << nsx*nsy << ", Total: " << nsx*nsy*ndx*ndy << endl;
    cout << "Run Time: " << runTime << " (" << printAsTime(runTime) << "), Sim Time: " << time << endl;
    cout << "Start Time: " << simulator.getStartRec() << ", Record Time: " << time - simulator.getStartRec() << endl;
    cout << "Iterations: " << iters << ", time per iter: " << (iters>0 ? toStr(runTime/iters) : "---") << endl;
    cout << "Transfer Time: " << transferTime << " (" << (runTime>0 ? toStr(transferTime/runTime*100) : "---") << "%)" << endl;
    cout << "Ratio: " << (runTime>0 ? toStr(time/runTime) : "---") << ", Ratio x Particles: " << (runTime>0 ? toStr(time/runTime*number) : "---") << endl;
    cout << " --- STATS ---\n";
    cout << "Neighbor list size: " << simulator.getNeighborListSize() << ", Ave per list: " << simulator.getAvePerNeighborList() << endl;
    cout << "----------------------- END SUMMARY -----------------------\n\n"; 

    /// Print recorded data
    if (animate) cout << simulator.printAnimationCommand(mode, novid) << endl;
    if (special) cout << simulator.printSpecialAnimationCommand(novid) << endl;
    if (bubbles) {
      cout << "bsize=" << simulator.getBubbleRecord() << ";\n"; //**
    }
    if (visBubbles) {
      if (!bubbles) cout << "bsize=" << simulator.getBubbleRecord() << ";\n";
      cout << "bubbles=" << simulator.getBubbles() << ";\n"; //**
    }
    string stats = simulator.printStatFunctions();
    if (!stats.empty()) cout << stats;
  }
  if (printSectors) simulator.printSectors();

  // End MPI
  MPI::Finalize();
  return 0;
}
