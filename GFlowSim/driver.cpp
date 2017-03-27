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
  bool animate = false;
  bool special = false;
  bool bubbles = false;
  bool omega   = false;
  bool KE      = false;
  bool LKE     = false;
  bool RKE     = false;
  bool cluster = false;
  bool triAlign = false;
  bool trackHeight  = false;
  bool novid   = false;
  bool printSectors = false;

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
  ArgParse parser(argc, argv);
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
  parser.get("animate", animate);
  parser.get("special", special);
  parser.get("bubbles", bubbles);
  parser.get("omega", omega);
  parser.get("KE", KE);
  parser.get("LKE", LKE);
  parser.get("RKE", RKE);
  parser.get("cluster", cluster);
  parser.get("triAlign", triAlign);
  parser.get("trackHeight", trackHeight);
  parser.get("novid", novid);
  parser.get("printSectors", printSectors);
  parser.get("square", square);
  parser.get("buoyancy", buoyancy);
  parser.get("loadFile", loadFile);
  parser.get("loadBuoyancy", loadBuoyancy);
  parser.get("saveFile", saveFile);
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
  if (0<skinDepth) simulator.setSkinDepth(skinDepth);
  // Create scenario
  if (loadFile=="" && loadBuoyancy=="") {
    if (buoyancy) simulator.createBuoyancyBox(radius, bR, density, width, height, velocity, dispersion);
    else if (square) simulator.createSquare(number, radius, width, height, velocity, dispersion);
    else throw false; // No selection
  }
  else if (loadBuoyancy!="") {
    
  }
  else { // We must have that loadFile!=""
    simulator.loadConfigurationFromFile(loadFile);
  }

  simulator.setTemperature(temperature);
  simulator.setDoInteractions(interact);
  simulator.setStartRec(start);

  simulator.setRecPositions(animate);
  simulator.setRecSpecial(special);
  simulator.setRecBubbles(bubbles);

  if (omega) simulator.addStatFunction(Stat_Omega, "omega");
  if (KE) simulator.addStatFunction(Stat_KE, "ke");
  if (LKE) simulator.addStatFunction(Stat_L_KE, "lke");
  if (RKE) simulator.addStatFunction(Stat_R_KE, "rke");
  if (cluster) simulator.addStatFunction(Stat_Clustering, "cluster");
  if (triAlign) simulator.addStatFunction(Stat_Triangle_Align, "triAlign");
  if (trackHeight) simulator.addStatFunction(Stat_Large_Object_Height, "height");

  // Get the actual number of particles in the simulation
  number = simulator.getSize();
  if (rank==0) {
    cout << "Number: " << number << endl;
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
    cout << "----------------------- END SUMMARY -----------------------\n\n"; 

    /// Print recorded data
    if (animate) cout << simulator.printAnimationCommand(novid) << endl;
    if (special) cout << simulator.printSpecialAnimationCommand(novid) << endl;
    if (bubbles) {
      cout << "bsize=" << simulator.getBubbleRecord() << ";\n"; //**
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
