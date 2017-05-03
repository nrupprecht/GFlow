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
  double density = 100;
  double velocity = 0.25;
  double omega = 0;
  double dispersion = 0;
  double temperature = 0;
  double gravity = -1e9;
  double drag = 0;
  double coeff = -1;
  double time = 1.;
  double start = 0;
  double epsilon = 1.e-4;
  double phi = 0.5;
  double skinDepth = -1;
  int interaction = 0;
  int LJ = -1;
  int Tri = -1;
  bool interact = true;
  bool seedRand = true;
  int quiet     = -1; // Quiet==-1 -> Normal, Quiet==1 -> Full quiet, Quiet==2 -> Quiet if label!=0
  int lattice   = 0;

  // Animation Paramaters
  int mode     = 0;     // Animation mode
  bool noprint = false;
  bool animate = false;
  int center   = false;
  bool snapshot = false;
  bool special = false;
  bool forces  = false;
  bool pressure = false;
  int forceChoice = 0;
  int typeChoice = 0;
  bool bubbles = false;
  bool visBubbles = false;
  bool bubbleField = false;
  bool csv = true;
  bool writeFields = false;
  bool writeFitness = false;
  bool writeAnimation = true;
  bool writeCreation = false;
  bool bulk    = false;
  // Stat functions
  bool angular = false;
  bool KE      = false;
  bool KEX     = false;
  bool KEY     = false;
  bool LKE     = false;
  bool RKE     = false;
  bool cluster = false;
  bool triAlign = false;
  bool trackHeight  = false;
  bool trackX  = false;
  bool GPE     = false;
  bool maxV    = false;
  bool maxVY   = false;
  bool largeVX = false;
  bool largeVY = false;
  bool num     = false;
  // Stat plots
  bool velDist = false;
  bool pressurePlot = false;
  bool densityPlot = false;
  // Options
  bool novid   = false;
  bool printSectors = false;
  double fps = -1; // FPS
  string label = "";
  double scale = 100;

  // Simulation type
  bool square = true;
  bool buoyancy = false;
  bool CV = false;
  bool CO = false;
  bool bacteria = false;

  // Load and save
  string loadFile = "";
  string loadBuoyancy = "";
  string createTube = "";
  string loadTube = "";
  string saveFile = "";
  string writeDirectory = "RunData";
  string updateBuoyancy = "";

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
  parser.get("omega", omega);
  parser.get("dispersion", dispersion);
  parser.get("temperature", temperature);
  parser.get("gravity", gravity);
  parser.get("drag", drag);
  parser.get("coeff", coeff);
  parser.get("time", time);
  parser.get("start", start);
  parser.get("epsilon", epsilon);
  parser.get("phi", phi);
  parser.get("skinDepth", skinDepth);
  parser.get("interaction", interaction);
  parser.get("LJ", LJ);
  parser.get("Tri", Tri);
  parser.get("interact", interact);
  parser.get("srand", seedRand);
  parser.get("quiet", quiet);
  parser.get("lattice", lattice);
  parser.get("mode", mode);
  parser.get("noprint", noprint);
  parser.get("animate", animate);
  parser.get("center", center);
  parser.get("snapshot", snapshot);
  parser.get("special", special);
  parser.get("forces", forces);
  parser.get("pressure", pressure);
  parser.get("forceChoice", forceChoice);
  parser.get("typeChoice", typeChoice);
  parser.get("bubbles", bubbles);
  parser.get("visBubbles", visBubbles);
  parser.get("bubbleField", bubbleField);
  parser.get("csv", csv);
  parser.get("writeFields", writeFields);
  parser.get("writeFitness", writeFitness);
  parser.get("writeAnimation", writeAnimation);
  parser.get("writeCreation", writeCreation);
  parser.get("bulk", bulk);
  parser.get("angular", angular);
  parser.get("KE", KE);
  parser.get("KEX", KEX);
  parser.get("KEY", KEY);
  parser.get("LKE", LKE);
  parser.get("RKE", RKE);
  parser.get("cluster", cluster);
  parser.get("triAlign", triAlign);
  parser.get("trackHeight", trackHeight);
  parser.get("trackX", trackX);
  parser.get("GPE", GPE);
  parser.get("maxV", maxV);
  parser.get("maxVY", maxVY);
  parser.get("largeVX", largeVX);
  parser.get("largeVY", largeVY);
  parser.get("num", num);
  parser.get("velDist", velDist);
  parser.get("pressurePlot", pressurePlot);
  parser.get("densityPlot", densityPlot);
  parser.get("novid", novid);
  parser.get("printSectors", printSectors);
  parser.get("fps", fps);
  parser.get("label", label);
  parser.get("scale", scale);
  parser.get("square", square);
  parser.get("buoyancy", buoyancy);
  parser.get("CV", CV);
  parser.get("CO", CO);
  parser.get("bacteria", bacteria);
  parser.get("loadFile", loadFile);
  parser.get("loadBuoyancy", loadBuoyancy);
  parser.get("createTube", createTube);
  parser.get("loadTube", loadTube);
  parser.get("saveFile", saveFile);
  parser.get("writeDirectory", writeDirectory);
  parser.get("updateBuoyancy", updateBuoyancy);
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
  if ((rank==0 && quiet==-1) || (quiet==2 && label=="0")) {
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
  simulator.setWriteCreation(writeCreation);
  simulator.setEpsilon(epsilon);
  simulator.setLatticeType(lattice);
  if (0<fps) simulator.setFrameRate(fps);
  if (0<skinDepth) simulator.setSkinDepth(skinDepth);
  // Adjust interaction
  if (LJ!=-1) interaction = 1;
  else if (Tri!=-1) interaction = 2;
  // Create scenario
  if (loadFile=="" && loadBuoyancy=="" && updateBuoyancy=="" && createTube=="" && loadTube=="") {
    if (buoyancy) simulator.createBuoyancyBox(radius, width, height, dispersion, interaction);
    else if (square) simulator.createSquare(number, radius, width, height, velocity, dispersion, interaction);
    else throw false; // No selection
  }
  else if (updateBuoyancy!="") {
    // Use bR = 0
    if (!simulator.loadBuoyancy(updateBuoyancy, 0, velocity, density)) { 
      cout << "Failed to load [" << updateBuoyancy << "], exiting.\n";
      return 0;
    }
    if (LJ!=-1) simulator.setInteractionType(1);
  }
  else if (loadBuoyancy!="") {
    if (!simulator.loadBuoyancy(loadBuoyancy, bR, velocity, density, CV, omega, CO)) {
      cout << "Failed to load [" << loadBuoyancy << "], exiting.\n";
      return 0;
    }
    if (LJ!=-1) simulator.setInteractionType(1);
  }
  else if (createTube!="") {
    if (!simulator.createTube(createTube)) {
      cout << "Failed to load [" << createTube << "], exiting.\n";
      return 0;
    }
  }
  else if (loadTube!="") {
    if (!simulator.loadTube(loadTube, bR, velocity, density, CV, omega, CO)) {
      cout << "Failed to load [" << loadTube << "], exitiny.\n";
      return 0;
    }
  }
  else { // We must have that loadFile!=""
    if (!simulator.loadConfigurationFromFile(loadFile)) {
      cout << "Failed to load [" << loadFile <<"], exiting.\n";
      return 0;
    }
    if (LJ!=-1) simulator.setInteractionType(1);
  }

  if (bacteria) simulator.setAsBacteria();

  simulator.setTemperature(temperature);
  if (gravity!=-1e9) simulator.setGravity(vec2(0,gravity));
  simulator.setDrag(drag);
  if (0<=coeff) simulator.setCoeff(coeff);
  simulator.setDoInteractions(interact);
  simulator.setStartRec(start);

  if (animate) simulator.setPositionsOption(1);
  if (center)  simulator.setFollowBall(true);
  simulator.setRecSpecial(special);
  if (forces)  simulator.setPositionsOption(2);
  simulator.setRecPressure(pressure);
  simulator.setForceChoice(forceChoice);
  simulator.setTypeChoice(typeChoice);
  simulator.setRecBubbles(bubbles);
  simulator.setVisBubbles(visBubbles);
  simulator.setBubbleField(bubbleField);
  simulator.setCSV(csv);
  // Directories
  simulator.setWriteFields(writeFields);
  simulator.setWriteFitness(writeFitness);
  simulator.setWriteAnimation(writeAnimation);
  simulator.setWriteDirectory(writeDirectory);
  // simulator.setRecBulk(bulk);
  if (buoyancy || loadBuoyancy!="") simulator.setFollowBall(true);
  // Stat functions
  if (angular) simulator.addStatFunction(Stat_Omega, "omega");
  if (KE) simulator.addStatFunction(Stat_KE, "ke");
  if (KEX) simulator.addStatFunction(Stat_KE_X, "kex");
  if (KEY) simulator.addStatFunction(Stat_KE_Y, "key");
  if (LKE) simulator.addStatFunction(Stat_L_KE, "lke");
  if (RKE) simulator.addStatFunction(Stat_R_KE, "rke");
  if (cluster) simulator.addStatFunction(Stat_Clustering, "cluster");
  if (triAlign) simulator.addStatFunction(Stat_Triangle_Align, "triAlign");
  if (trackHeight) simulator.addStatFunction(Stat_Large_Object_Height, "height");
  if (trackX) simulator.addStatFunction(Stat_Large_Object_X, "posx");
  if (GPE) simulator.addStatFunction(Stat_Gravitational_PE, "gpe");
  if (maxV) simulator.addStatFunction(Stat_Max_Velocity, "maxv");
  if (maxVY) simulator.addStatFunction(Stat_Max_Velocity_Y, "maxvy");
  if (largeVX) simulator.addStatFunction(Stat_Large_Object_VX, "largevx");
  if (largeVY) simulator.addStatFunction(Stat_Large_Object_VY, "largevy");
  if (num) simulator.addStatFunction(Stat_Number_Particles, "num");
  // Stat plots
  if (velDist) simulator.addStatPlot(Plot_Velocity, "velDist", 0, 3); //** 3 is a magic number
  if (pressurePlot) simulator.addStatPlot(Plot_Force_Vs_Depth, "pressPlot", simulator.getBottom(), simulator.getTop());
  if (densityPlot) simulator.addStatPlot(Plot_Particle_Density_Vs_Depth, "denPlot", simulator.getBottom(), simulator.getTop());

  // Get the actual number of particles in the simulation
  number = simulator.getSize();
  double filled = simulator.getFilledVolume(), vol = simulator.getVolume();
  if ((rank==0 && quiet==-1) || (quiet==2 && label=="0")) {
    cout << "Number: " << number << ", Packing fraction: " << filled/vol << endl;
    cout << "Dimensions: " << simulator.getWidth() << " x " << simulator.getHeight() << endl;
    cout << "Set up time: " << simulator.getSetUpTime() << endl;
    cout << "  ..........................\n";
  }

  if (interaction!=0) simulator.setInteractionType(interaction);

  // Run the simulation
  simulator.run(time);

  // Head node prints the run summary
  if ((rank==0 && quiet==-1) || (quiet==2 && label=="0")) {
    if (updateBuoyancy!="") {
      if (simulator.createConfigurationFile(updateBuoyancy)) cout << "Saved configuration to file [" << updateBuoyancy << "]" << endl;
    }
    if (saveFile!="") {
      if (simulator.createConfigurationFile(saveFile)) cout << "Saved configuration to file [" << saveFile << "]" << endl;
    }

    double runTime = simulator.getRunTime(), transferTime = simulator.getTransferTime();
    int iters = simulator.getIter(), ndx = simulator.getNDX(), ndy = simulator.getNDY();
    int nsx = simulator.getNSX(), nsy = simulator.getNSY();
    cout << "Domains: " << ndx << " x " << ndy << ", Total: " << ndx*ndy << endl;
    cout << "Sectors: " << nsx << " x " << nsy << ", Per Domain: " << nsx*nsy << ", Total: " << nsx*nsy*ndx*ndy << endl;
    cout << "Run Time: " << runTime << " (" << printAsTime(runTime) << "), Sim Time: " << simulator.getTime() << endl;
    cout << "Start Time: " << simulator.getStartRec() << ", Record Time: " << time - simulator.getStartRec() << endl;
    cout << "Iterations: " << iters << ", time per iter: " << (iters>0 ? toStr(runTime/iters) : "---") << endl;
    double transferPercent = transferTime/runTime*100;
    cout << "Transfer Time: " << transferTime << " (" << (runTime>0 ? (transferPercent>0.01 ? toStr(transferPercent) : "~ 0") : "---") << "%)" << endl;
    cout << "Ratio: " << (runTime>0 ? toStr(simulator.getTime()/runTime) : "---") << " (" << (runTime>0 ? toStr(runTime/time) : "---") <<  "), Ratio x Particles: " << (runTime>0 ? toStr(time/runTime*number) : "---") << endl;
    cout << " --- STATS ---\n";
    cout << "Neighbor list size: " << simulator.getNeighborListSize() << ", Ave per list: " << simulator.getAvePerNeighborList() << endl;
    cout << "----------------------- END SUMMARY -----------------------\n\n"; 
  }
  if (rank==0) { // Do this even when quiet = true
    /// Print recorded data
    simulator.setScale(scale);
    if (animate && !writeAnimation) cout << simulator.printAnimationCommand(mode, novid, label) << endl;
    if (snapshot) cout << simulator.printSnapshot() << endl;
    if (special) cout << simulator.printSpecialAnimationCommand(novid) << endl;
    if (forces && !writeAnimation)  cout << simulator.printForcesAnimationCommand(mode, novid) << endl;
    if (pressure) cout << simulator.getPressureRecord() << ";\n";
    // Print bubble data
    if (bubbles && !csv) {
      cout << "bsize" << label << "=" << simulator.getBubbleRecord() << ";\n";
      cout << "num" << label << "=Table[Length[bsize" << label << "[[i]]],{i,1,Length[bsize" << label << "]}];\n";
      cout << "vol" << label << "=Table[Total[bsize" << label << "[[i]]],{i,1,Length[bsize" << label << "]}];\n";
      if (!noprint) {
        cout << "Print [\"Number of bubbles\"]\n";
        cout << "ListLinePlot[num" << label << ",PlotStyle->Black,ImageSize->Large,PlotRange -> All]\n";
        cout << "Print[\"Total bubble volume\"]\n";
        cout << "ListLinePlot[vol" << label << ",PlotStyle->Black,ImageSize->Large,PlotRange->All]\n";
      }
    }
    if (visBubbles) cout << simulator.printBulkAnimationCommand(novid) << endl;
    if (bubbleField && !csv) cout << "bubF=" << mmPreproc(simulator.getBubbleField()) << ";\n";
    // Print statistics and plots
    string stats = simulator.printStatFunctions(label, noprint);
    if (!stats.empty()) cout << stats;
    string plots = simulator.printStatPlots(label, noprint);
    if (!plots.empty()) cout << plots;
  }
  if (printSectors) simulator.printSectors();

  // End MPI
  MPI::Finalize();
  return 0;
}
