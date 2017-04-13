#ifndef GFLOWBASE_H
#define GFLOWBASE_H

#include "Sectorization.h"
#include <functional>
#include <tuple> 

typedef pair<vec2, vec2> vpair;
//typedef std::tuple<vec2, double, double> Tri;

// Zero velocity and angular velocity functions
inline vec2 ZeroV(double) { return 0; }
inline double ZeroOm(double) { return 0; }

// Granular float base class
class GFlowBase {
 public:
  GFlowBase();
  ~GFlowBase();

  void run(double);  // Run the simulation for some amount of time

  // Addition
  void addWall(Wall);
  void addWall(double, double, double, double);
  void addParticle(Particle);
  void addParticle(double, double, double);

  // Accessors
  Bounds getBounds()       { return Bounds(left, right, bottom, top); }
  int getNSX()             { return sectorization.getNSX(); }
  int getNSY()             { return sectorization.getNSY(); }
  int getNDX()             { return ndx; }
  int getNDY()             { return ndy; }
  int getSize();
  double getFilledVolume();
  double getVolume()       { return (right-left)*(top-bottom); }
  double getWidth()        { return right-left; }
  double getHeight()       { return top-bottom; }
  bool getWrapX()          { return wrapX; }
  bool getWrapY()          { return wrapY; }
  vec2 getGravity()        { return gravity; }
  double getTemperature()  { return temperature; }
  double getViscosity()    { return viscosity; }
  double getTime()         { return time; }
  double getEpsilon()      { return epsilon; }
  double getDispTime()     { return dispTime; }
  double getDispRate()     { return 1./dispTime; }
  double getStartRec()     { return startRec; }
  int getRecIter()         { return recIter; }
  int getIter()            { return iter; }
  double getRunTime()      { return runTime; }
  bool getRunning()        { return running; }
  double getTransferTime() { return transferTime+sectorization.getTransferTime(); }
  double getCutoff()       { return sectorization.getCutoff(); }
  bool getDoInteractions() { return doInteractions; }

  // Statistics accessors
  int getNeighborListSize() { return sectorization.getNeighborListSize(); }
  double getAvePerNeighborList() { return sectorization.getAvePerNeighborList(); }

  // Mutators
  void setBounds(double, double, double, double);
  void setBounds(Bounds);
  void setWrapX(bool w)          { wrapX = w; }
  void setWrapY(bool w)          { wrapY = w; }
  void setGravity(vec2 g);
  void setTemperature(double t);
  void setDrag(double d)         { sectorization.setDrag(d); }
  void setViscosity(double h);
  void setLatticeType(int l)     { latticeType = l; }
  void setStartRec(double s)     { startRec = s; }
  void setDoInteractions(bool i);
  void setInteractionType(int i);
  void setExpectedSize(int i)    { sectorization.setASize(i); }
  void setSkinDepth(double s)    { skinDepth = s; sectorization.setSkinDepth(s); }
  void setFrameRate(double f)    { dispTime = 1./f; }
  void setEpsilon(double);

  // File functions
  virtual bool loadConfigurationFromFile (string);
  virtual bool createConfigurationFile   (string);

  // -----  TO GO TO GFLOW.H  ------
  void setAsBacteria();

  void createSquare(int, double, double=4., double=4., double=0.1, double=0., int=0);
  void createBuoyancyBox(double,double,double,double,double,double,double, int=0);
  bool loadBuoyancy(string, double=0.5, double=5., double=10.);
  void recordPositions();

  void addStatFunction(StatFunc, string);
  string printStatFunctions(string="");

  string printWallsCommand();
  string printAnimationCommand(int=0, bool=false, string="");
  string printSpecialAnimationCommand(bool=false);
  string printForcesAnimationCommand(int=0, bool=false);
  string printBulkAnimationCommand(bool=false);
  void printSectors();

  // Get bubble sizes variants and helper functions
  void getBulkData(vector<VPair>&, double=0.01, double=0.025, double=0.05);
  vector<double> getBulkData(vector<Particle>&, Bounds=NullBounds, double=0.01, double=0.0001, double=0.01);
  vector<double> getBulkData(vector<Particle>&, vector<VPair>&, Bounds=NullBounds, double=0.01, double=0.0001, double=0.01);
  vector<double> getBulkData(vector<Particle>&, string&, vector<VPair>&, bool, Bounds=NullBounds, double=0.01, double=0.0001, double=0.01);
  inline void unite(int*, int, int);
  inline int getHead(int*, int);
  inline void createOutline(int*, int, int, double, double, Bounds, vector<VPair>&);
  inline void createMatrix(int*, int, int, double, double, double, vector<int>, string&);

  // Record accessors
  auto getPositionRecord() { return positionRecord; }
  auto getSpecialRecord()  { return specialRecord; }
  auto getForceRecord()    { return forceRecord; }
  auto getBubbleRecord()   { return bubbleRecord; }
  auto getBulkRecord()     { return bulkRecord; }
  string printPositionRecord(int);
  vector<vpair> getWallsPositions();

  double getSetUpTime() { return setUpTime; }

  void setRecPositions(bool b) { recPositions = b; }
  void setRecSpecial(bool b)   { recSpecial = b; }
  void setRecForces(bool b)    { recForces = b; }
  void setRecBubbles(bool b)   { recBubbles = b; }
  void setVisBubbles(bool b)   { visBubbles = b; }
  void setRecBulk(bool b)      { recBulk = b; }
  void setRestrictBubbleDomain(bool b) { restrictBubbleDomain = b; }
  void setForceChoice(int c)   { forceChoice = c; }
 protected:
  double setUpTime;
  // Data
  vector<vector<PData> > positionRecord;
  vector<vector<Tri> > specialRecord;
  vector<vector<Tri> > forceRecord;
  vector<vector<double> > bubbleRecord;
  vector<vector<VPair> > bulkRecord;
  bool recPositions;
  bool recSpecial;
  bool recForces;
  bool recBubbles;
  bool visBubbles;
  bool recBulk;
  bool restrictBubbleDomain;
  int forceChoice;

  vector<pair<StatFunc,string> > statFunctions; // Statistic functions and a string to name them
  vector<vector<vec2> >  statRecord;    // Save the data produced by the statistic functions

 public:
  // -------------------------------

 protected:
  /// Principal functions
  virtual void setUpSectorization(); // Set up the sectorization
  virtual void setUpSectorization(Sectorization&, vec2); // Set up a sectorization
  virtual void resetVariables();     // Reset variables for the start of a simulation
  virtual void objectUpdates();      // Do forces, move objects
  virtual void logisticUpdates();    // Update times
  virtual void record();             // Record data
  virtual void checkTermination();   // Check whether we should terminate the program
  virtual void resets();             // Reset objects as neccessary
  virtual void gatherData();         // Gather data back to processor 0
  virtual void discard();            // Discard and reset the simulation state

  /// Helper functions
  void bestProcessorGrid(int&, int&, const int, const Bounds);
  Bounds getBoundsForProc(int);
  Bounds getBoundsForProc(int, const Bounds&);
  void distributeParticles(list<Particle>&, Sectorization&);
  void recallParticles(vector<Particle>&); // Gather copies of all particles into a vector
  void recallParticlesByProcessor(vector<vector<Particle> >&); // Gather copies of all particles, remembering which particles came from which processor
  double reduceStatFunction(StatFunc, int=0);

  list<Particle> createParticles(vector<vec2>, double, double, std::function<vec2(double)> = ZeroV, std::function<double(double)> = ZeroOm, double=default_sphere_coeff, double=default_sphere_dissipation, double=default_sphere_repulsion, int=0);
  // Create and distribute particles variants
  void createAndDistributeParticles(int, const Bounds&, Sectorization&, double, double=0, std::function<vec2(double)> = ZeroV, double=default_sphere_coeff, double=default_sphere_dissipation, double=default_sphere_repulsion, int=0);
  void createAndDistributeParticles(const vector<vec2>&, const vector<double>&, const vector<int>&, const Bounds&, Sectorization&, std::function<vec2(double)> = ZeroV, double=default_sphere_coeff, double=default_sphere_dissipation, double=default_sphere_repulsion);
  void createAndDistributeParticles(const vector<double>&, const vector<int>&, const Bounds&, Sectorization&, std::function<vec2(double)> = ZeroV, double=default_sphere_coeff, double=default_sphere_dissipation, double=default_sphere_repulsion);
  // Find packed solution variants
  vector<vec2> findPackedSolution(int, double, Bounds, vec2 = 0, double=0.5, double=0.5);
  vector<vec2> findPackedSolution(const vector<double>&, const vector<int>&, const Bounds&, vec2=0, double=0.5, double=0.5);
  // Find lattice solution
  vector<vec2> findLatticeSolution(int, double, Bounds, int=-1, double=0.1, double=0.2);

  /// Data
  double left, right, bottom, top;
  bool wrapX, wrapY;
  vec2 gravity;
  double temperature, viscosity;

  /// Times
  double time;                   // Simulated time
  double epsilon;                // Time step and its square root
  double dispTime, lastDisp;     // Time between recording (1. / display rate), last time data was recorded
  double startRec;               // When to start recording data
  int iter, recIter, maxIter;    // Simulation iteration, how many iterations we have recorded data at, maximum iteration
  double runTime;                // How much (real) time the last simulation took to run
  bool running;                  // True if the simulation is currently runnint
  double transferTime;           // How much time is spent by MPI transfering data
  // Initialization
  int latticeType;               // What type of lattice we should initialize, -1 means use find packed solution

  /// Objects
  vector<Wall> walls;            // A vector of all the walls in the simulation
  list<Particle> particles;      // A vector of all the particles in the simulation
  bool doInteractions;           // True if we let the particles interact with one another

  /// Termination
  vector<pair<StatFunc, std::function<bool(double)> > > terminationConditions;
  double startChecking, lastCheck, checkDelay;

  /// Sectorization
  Sectorization sectorization;   // The sectorization for this processor
  double skinDepth;           // What skin depth we should use for our sectors
  bool doWork;                   // True if this processor needs to do work

  // MPI
  int rank, numProc;             // The rank of this processor and the total number of processors
  int ndx, ndy;                  // Number of domains we divide into
  MPI_Datatype PARTICLE;
  MPI::Intercomm CommWork;       // The communicator for the working processors
};

#endif
