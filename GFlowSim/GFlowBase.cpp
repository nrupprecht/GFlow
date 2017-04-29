#include "GFlowBase.h"

using MPI::COMM_WORLD;

GFlowBase::GFlowBase() {
  left = 0; right = 1; bottom = 0; top = 1;
  wrapX = false; wrapY = false;
  gravity = vec2(0., -1.);
  temperature = 0; viscosity = 1.308e-3;
  latticeType = 0;
  time = 0;
  epsilon = 1e-4;
  sectorization.setEpsilon(epsilon);
  dispTime = 1./15.; lastDisp = 0;
  startRec = 0;
  iter = 0; recIter = 0; maxIter = 0;
  runTime = 0;
  running = false;
  transferTime = 0;
  doInteractions = true;
  startChecking = 1.;
  checkDelay = 1./15.; lastCheck = 0.;
  skinDepth = 0.025;
  doWork = false;

  //---
  setUpTime = 0;
  recSpecial = false;
  visBounds = bubbleBounds = NullBounds;
  writeFields = false;
  writeFitness = false;
  writeAnimation = true; // We animate by printing files by default
  writeCreation = false;
  forceChoice = 0;
  typeChoice = 0;
  doFields = false;
  fieldUpdateDelay = 0.0005;
  fieldUpdateCounter = 0;
  scale = 100;
  writeDirectory = "RunData";
  statPlotBins = 100;
  alphaR = default_alphaR;
  alphaW = default_alphaW;
  csatR = default_csatR;
  csatW = default_csatW;
  betaR = default_betaR;
  for (int i=0; i<=7; ++i) options[i] = 0;
  //---

  // Get MPI system data
  rank = COMM_WORLD.Get_rank();
  numProc = COMM_WORLD.Get_size();
  CommWork = COMM_WORLD;

  // Define the particle data type
  MPI_Type_contiguous( 16*sizeof(double)/sizeof(float), MPI_FLOAT, &PARTICLE );
  MPI_Type_commit( &PARTICLE );
}

GFlowBase::~GFlowBase() {}

void GFlowBase::run(double runLength) {
  setUp();
  resetVariables();
  // Create work communicator
  //** int color = doWork ? 1 : 0;
  //** MPI::Intercomm CommWork = COMM_WORLD.Split(color, rank);
  sectorization.setCommWork(CommWork);
  // Calculate the number of iterations we will run for
  maxIter = runLength/epsilon;
  auto start = current_time();
  running = true;
  if (startRec<=0) record(); // Initial record
  if (doWork) {
    for (iter=0; iter<maxIter && running; ++iter) {
      objectUpdates();
      if (time-lastDisp>dispTime && startRec<time) record();
      logisticUpdates();
      if (time-lastCheck>checkDelay && time>startChecking) checkTermination();
      CommWork.Barrier();
    }
  }
  running = false;
  auto end = current_time();
  runTime = time_span(end, start);
  // Release work comm and set sectorization comm back to comm world
  sectorization.resetComm();
  // CommWork.Free(); 
  // Update transfer time
  transferTime += sectorization.getTransferTime();
  endOfRun();
}

void GFlowBase::addWall(Wall w) {
  walls.push_back(w);
  sectorization.addWall(w);
}

void GFlowBase::addWall(double sx, double sy, double ex, double ey) {
  Wall w(vec2(sx,sy), vec2(ex,ey));
  walls.push_back(w);
  sectorization.addWall(w);
}

void GFlowBase::addParticle(Particle p) {
  particles.push_back(p);
}

void GFlowBase::addParticle(double x, double y, double r) {
  particles.push_back(Particle(x, y, r));
}

int GFlowBase::getSize() {
  int sz = sectorization.getSize();
  int size = 0;
  CommWork.Reduce(&sz, &size, 1, MPI_INT, MPI_SUM, 0);
  return size;
}

double GFlowBase::getFilledVolume() {
  return reduceStatFunction(Stat_Volume, 2);
}

void GFlowBase::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
  sectorization.setSimBounds(l, r, b, t);
}

void GFlowBase::setBounds(Bounds b) {
  left = b.left; right = b.right; bottom = b.bottom; top = b.top;
  sectorization.setSimBounds(b);
}

void GFlowBase::setGravity(vec2 g) {
  gravity = g;
  sectorization.setGravity(g);
}

void GFlowBase::setTemperature(double t) {
  temperature = t;
  sectorization.setTemperature(t);
}

void GFlowBase::setViscosity(double h) {
  viscosity = h;
  sectorization.setViscosity(h);
}

void GFlowBase::setDoInteractions(bool i) {
  doInteractions = i;
  sectorization.setDoInteractions(i);
}

void GFlowBase::setInteractionType(int i) {
  sectorization.setInteractionType(i);
}

void GFlowBase::setEpsilon(double ep) {
  epsilon = ep;
  sectorization.setEpsilon(ep);
}

bool GFlowBase::loadConfigurationFromFile (string filename) {
  discard();
  if (rank==0) {
    vect<> g;
    vector<double> radii;
    vector<vect<> > wallLeft, wallRight, positions;
    vector<int> interactions;
    ifstream fin(filename);
    if (fin.fail()) return false;
    // Get simulation bounds
    fin >> left >> right >> bottom >> top;
    // Get gravity
    fin >> g;
    // Get walls
    fin >> wallLeft;
    fin >> wallRight;
    // Get radii and positions
    fin >> radii;
    fin >> positions;
    fin >> interactions;
    fin.close();
    // Clear current configuration
    discard();
    // Set up
    setBounds(left, right, bottom, top);
    gravity = g;
    int size = positions.size(), rsize = radii.size(), wsize = min(wallLeft.size(), wallRight.size()), isize = interactions.size();
    for (int i=0; i<wsize; ++i)
      addWall(Wall(wallLeft[i], wallRight[i]));
    for (int i=0; i<size; ++i) {
      Particle p(positions.at(i), radii.at(i%rsize));
      p.interaction = interactions.at(i%isize);
      addParticle(p);
    }
    // Send out particles
    setUpSectorization();
    distributeParticles(particles, sectorization);
    sectorization.initialize();
    // Return success
    return true;
  }
  return false;
}

bool GFlowBase::createConfigurationFile (string filename) {
  if (rank==0) {
    vector<double> radii;
    vector<vec2> wallLeft, wallRight, positions;
    vector<int> interactions;
    // Accumulate data
    ofstream fout(filename);
    if (fout.fail()) return false;
    // Recall particles
    vector<Particle> allParticles;
    recallParticles(allParticles);
    // Record wall positions
    for (const auto &W : walls) {
      wallLeft.push_back(W.getLeft());
      wallRight.push_back(W.getRight());
    }
    // Record particle positions
    double sig = 0;
    bool sameSig = true, sameInteraction = true;
    int inter = 0;
    if (!allParticles.empty()) {
      sig = allParticles.at(0).sigma;
      inter = allParticles.at(0).interaction;
      for (const auto &P : allParticles) {
	radii.push_back(P.sigma);
	if (P.sigma!=sig) sameSig = false;
	positions.push_back(P.position);
	interactions.push_back(P.interaction);
	if (P.interaction!=inter) sameInteraction = false;
      }
    }
    // Write data to file
    fout << left << " " << right << " " << bottom << " " << top << "\n";
    fout << gravity << "\n";
    // Print wall data
    fout << wallLeft << "\n";
    fout << wallRight << "\n";
    // Print radii (or radius if they are all the same)
    if (sameSig) fout << "{" << sig << "}\n";
    else fout << radii << "\n";
    // Print positions
    fout << positions << "\n";
    // Print interactions
    if (sameInteraction) fout << "{" << inter << "}\n";
    else fout << interactions << "\n";
    fout.close();
    return true;
  }
  return false;
}

void GFlowBase::setUpSectorization() {
  // Decide how to divide up the space into domains
  bestProcessorGrid(ndx, ndy, numProc, Bounds(left, right, bottom, top));
  sectorization.giveDomainInfo(ndx, ndy);
  // Calculate what bounds this processor is in charge of
  if (rank<ndx*ndy) {
    doWork = true; // This processor needs to do work
    sectorization.setSimBounds(left, right, bottom, top);
    sectorization.setSkinDepth(skinDepth);
    Bounds domainBounds = getBoundsForProc(rank);
    sectorization.setBounds(domainBounds);
    // sectorization.setASize(particles.size()); //** AD HOC
    sectorization.setGravity(gravity);
    sectorization.setDoInteractions(doInteractions);
  }
  else doWork = false;
}

void GFlowBase::setUpSectorization(Sectorization &sectors, vec2 grav) {
  // Decide how to divide up the space into domains
  bestProcessorGrid(ndx, ndy, numProc, sectors.getSimBounds() );
  sectors.giveDomainInfo(ndx, ndy);
  // Calculate what bounds this processor is in charge of
  if (rank<ndx*ndy) {
    doWork = true; // This processor needs to do work
    sectors.setSimBounds(left, right, bottom, top);
    sectors.setSkinDepth(skinDepth);
    Bounds domainBounds = getBoundsForProc(rank, sectors.getSimBounds() );
    sectors.setBounds(domainBounds);
    sectors.setGravity(grav);
    sectors.setDoInteractions(doInteractions);
    sectors.initialize(); // Keep this here (for now)
  }
  else doWork = false;
}

void GFlowBase::setUp() {
  if (writeFields || writeFitness || writeAnimation || writeCreation) {
    // Remove any previously existing file
    system(("rm -r "+writeDirectory).c_str());
    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);
    if (writeFields) {
      mkdir((writeDirectory+"/Waste").c_str(), 0777); // Waste director
      mkdir((writeDirectory+"/Resource").c_str(), 0777); // Resource directory
    }
    if (writeFitness) mkdir((writeDirectory+"/Fitness").c_str(), 0777); // Fitness director
    if (writeAnimation) mkdir((writeDirectory+"/Pos").c_str(), 0777); // Position directory
    if (writeCreation) mkdir((writeDirectory+"/Init").c_str(), 0777); // Initialization directory
  }
  resetVariables();
}

void GFlowBase::resetVariables () {
  time = 0;
  iter = 0;
  recIter = 0;
  maxIter = 0;
  runTime = 0;
}

void GFlowBase::objectUpdates() {
  // Update particles
  sectorization.update();
  // Update fields
  if (doFields) {
    // Get the particle
    double *px = sectorization.getPX();
    double *py = sectorization.getPY();
    int *it = sectorization.getIT();
    Characteristic **ch = sectorization.getCH();
    int array_end = sectorization.getArrayEnd();
    // Eat and produce waste
    for (int i=0; i<array_end; ++i) 
      if (-1<it[i]) {
	Bacteria *b = reinterpret_cast<Bacteria*>(ch[i]);
	if (b!=0) {
	  double rSec = b->secretion;	  
	  Resource.increase(px[i], py[i], rSec*epsilon);
	  Waste.increase(px[i], py[i], default_bacteria_waste*epsilon);
	}
      }
    // Update fields
    if (fieldUpdateDelay<fieldUpdateCounter) {
      Resource.update(epsilon);
      Waste.update(epsilon);
      fieldUpdateCounter = 0;
    }
    fieldUpdateCounter += epsilon;
    // Set fitnesses
    for (int i=0; i<array_end; ++i)
      if (-1<it[i]) {
	Bacteria *b = reinterpret_cast<Bacteria*>(ch[i]);
	if (b!=0) {
	  double res = Resource.at(px[i],py[i]), wst = Waste.at(px[i],py[i]);
	  double rSec = b->secretion;
	  double fitness = alphaR*res/(res+csatR) - alphaW*wst/(wst+csatW) - betaR*rSec;
	  //** There is a Heisenbug here. I don't know what it is.
	  b->setFitness(fitness);
	}
      }
  }
}

void GFlowBase::logisticUpdates() {
  time += epsilon;
}

void GFlowBase::record() {
  // So we only need to load these once
  vector<Particle> allParticles;
  // VISUALIZATION OPTIONS
  // [0] - Position animation options 1 -> Print positions, 2 -> Print pressures
  // [1] - Record number of bubbles
  // [2] - Record total bubble volume
  // [3] - Visualize bubbles (bulk animation)
  // [4] - Create bubble field
  // [5] - Record waste field
  // [6] - Record resource field
  // [7] - Record fitness field
  if (options[0]) {
    // Get the required data
    vector<PData> positions;
    if (options[0]==1) {
      recallParticles(allParticles);
      for (const auto& p : allParticles) positions.push_back(PData(p.position, p.sigma, p.theta, p.interaction, 0)); 
    }
    else { // Force or pressure data - ** Would have to do something different if using multiple processors
      positions = sectorization.forceAnimate(forceChoice, typeChoice);
    }
    // Write data
    if (writeAnimation)
      if (!printToCSV(writeDirectory+"/Pos/pos", positions, recIter))
	cout << "Printing to [" << writeDirectory << "/Pos/pos." << recIter << ".csv] Failed.\n";
    else positionRecord.push_back(positions);
  }
  // Bubble related options
  if (options[1] || options[2] || options[3] || options[4]) {
    // Find the proper bounds if neccessary
    if (followBall) bubbleBounds = followBallBounds();
    else if (bubbleBounds==NullBounds) bubbleBounds = Bounds(left, right, bottom, top);
    // Record what bounds we used
    animationBounds.push_back(bubbleBounds);
    // Get all bubble data
    vector<VPair> vis;
    if (allParticles.empty()) recallParticles(allParticles);
    auto bubbleData = getBulkData(allParticles, vis, bubbleBounds);
    bulkBounds.push_back(bubbleBounds);
    // Store required data
    if (options[1]) bubbleRecord.push_back(bubbleData);
    if (options[2]); // Created from bubbleRecord
    if (options[3]) bulkRecord.push_back(vis);
    if (options[4]); // Created automatically
  }
  // Print Waste field
  if (options[5] && !Waste.empty()) printToCSV(writeDirectory+"/Waste/wst", Waste, recIter);
  // Print Resource field
  if (options[6] && !Resource.empty()) printToCSV(writeDirectory+"/Resource/rsc", Resource, recIter);
  // Print Fitness field
  if (options[7] && !Waste.empty()) {
    updateFitness();
    printToCSV(writeDirectory+"/Fitness/fit", Fitness, recIter);
  }
  
  // Record stat function statistics
  if (!statFunctions.empty()) sectorization.updatePList();
  int i=0;
  for (auto &sf : statFunctions) {
    double data = reduceStatFunction(sf.first);
    statRecord.at(i).push_back(vec2(time, data));
    ++i;
  }
  // Record stat plot statistics
  i=0;
  for (auto &sp : statPlots) {
    sp.first(allParticles, statPlotRecord.at(i), statPlotBounds.at(i).first,  statPlotBounds.at(i).second);
  }
  // Special record
  if (recSpecial) {
    vector<vector<Particle> > specialView;
    recallParticlesByProcessor(specialView);
    int c = 0;
    vector<Tri> positions;
    for (const auto &lst : specialView) {
      for (const auto &p : lst)
	positions.push_back(Tri(p.position, p.sigma, c));
      ++c;
    } 
    specialRecord.push_back(positions);
  }
  // Update display 
  lastDisp = time;
  ++recIter;
}

void GFlowBase::checkTermination() {
  lastCheck = time;
  // Check termination conditions
  if (terminationConditions.empty()) return;
  vector<Particle> allParticles;
  recallParticles(allParticles);  
  for (const auto& pr : terminationConditions) {
    int count = 0;
    double data = pr.first(allParticles, count);
    if (pr.second(data)) {
      running = false;
      break;
    }
  }
}

void GFlowBase::resets() {

}

void GFlowBase::gatherData() {
  
}

void GFlowBase::endOfRun() {
  if ((writeAnimation || writeFields || writeFitness) && rank==0) {
    // Print a master file
    printToCSV(writeDirectory+"/number.csv", vector<int>(1,recIter)); // Print how many files to expect
    // Print walls to file
    printToCSV(writeDirectory+"/walls.csv", walls);
    // Print bounds to file
    Bounds bounds(left,right,bottom,top);
    printToCSV(writeDirectory+"/bnds.csv", vector<Bounds>(1,bounds));
  }
}

void GFlowBase::discard() {
  walls.clear();
  particles.clear();
  sectorization.discardAll();
}

void GFlowBase::bestProcessorGrid(int &x, int &y, const int number, const Bounds b) {
  int xm=1, ym=1;
  double rmin=1e9, r=0, l=(b.right-b.left)/(b.top-b.bottom);
  // Find the domain grid that is most square and uses the most number of processors
  for (int i=1; i<number; ++i) {
    int j=number/i;
    double factor = static_cast<double>(j)/static_cast<double>(i)*l;
    r = factor + 1./factor;
    if (r<rmin) { rmin = r; xm = i; ym = j; }
  } 
  x = xm; y = ym;
}

Bounds GFlowBase::getBoundsForProc(int rnk) {
  double dx = (right-left)/ndx, dy = (top-bottom)/ndy;
  double l = (rnk%ndx)*dx, r = l+dx;
  double b = (rnk/ndx)*dy, t = b+dy;
  return Bounds(l, r, b, t);
}

Bounds GFlowBase::getBoundsForProc(int rnk, const Bounds &bnds) {
  double dx = (bnds.right-bnds.left)/ndx, dy = (bnds.top-bnds.bottom)/ndy;
  double l = (rnk%ndx)*dx, r = l+dx;
  double b = (rnk/ndx)*dy, t = b+dy;
  return Bounds(l, r, b, t);
}

void GFlowBase::distributeParticles(list<Particle> &allParticles, Sectorization &sectors) {
  // Take your own particles
  list<Particle> parts;
  if (rank==0) {
    Bounds bnds = getBoundsForProc(0, sectors.getSimBounds() );
    vector<list<Particle>::iterator> remove;
    for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
      if (bnds.contains(p->position)) {
        remove.push_back(p);
        parts.push_back(*p);
      }
    for (auto &p : remove) allParticles.erase(p);
    // Add particles to sectors
    for (auto &p : parts) sectors.addParticle(p);
  }
  // Distribute particles to processes
  int max = min(ndx*ndy, numProc);
  for (int proc=1; proc<max; ++proc) {
    if (rank==0) { // Send
      vector<Particle> domainParticles;
      Bounds bnds = getBoundsForProc(proc, sectors.getSimBounds() );
      vector<list<Particle>::iterator> remove;
      for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
        if (bnds.contains(p->position)) {
          remove.push_back(p);
          domainParticles.push_back(*p);
        }
      for (auto &p : remove) allParticles.erase(p);
      int size = domainParticles.size();
      Particle *buffer = new Particle[size];
      for (int j=0; j<size; ++j) buffer[j] = domainParticles[j];
      // Send the amount of data we are going to send to processor i
      CommWork.Send( &size, 1, MPI_INT, proc, 0);
      // Send the actual data to processor i
      CommWork.Send( buffer, size, PARTICLE, proc, 0);
    }
    else if (rank==proc) { // Recieve
      int size = 0, root = 0;
      // Recieve the amount of data we should expect
      CommWork.Recv( &size, 1, MPI_INT, root, 0);
      // Recieve the actual data
      Particle *buffer = new Particle[size];
      CommWork.Recv( buffer, size, PARTICLE, root, 0);
      // Add particles to sectors
      for (int i=0; i<size; ++i) sectors.addParticle(buffer[i]);
    }
  }
  CommWork.Barrier();
}

void GFlowBase::recallParticles(vector<Particle>& allParticles) {
  // Get all the particles back from the sectorizations
  auto start = current_time();
  int max = min(ndx*ndy, numProc);
  if (rank==0)
    for (auto &p : sectorization.getParticles())
      allParticles.push_back(p);
  for (int proc=1; proc<max; ++proc) {
    if (rank==0) {
      int size = 0;
      // Recieve how many Particles to expect to recieve
      CommWork.Recv( &size, 1, MPI_INT, proc, 0);
      // Create a buffer of Particles to send
      Particle *buffer = new Particle[size];
      CommWork.Recv( buffer, size, PARTICLE, proc, 0);
      for (int i=0; i<size; ++i) allParticles.push_back(buffer[i]);
    }
    else if (rank==proc) {
      auto &parts = sectorization.getParticles();
      int size = parts.size(), root = 0;
      CommWork.Send( &size, 1, MPI_INT, root, 0);
      // Create a buffer of Particles to send
      Particle *buffer = new Particle[size];
      int i=0;
      for (auto p=parts.begin(); p!=parts.end(); ++p, ++i) buffer[i] = *p;
      CommWork.Send( buffer, size, PARTICLE, root, 0);
    }
  }
  CommWork.Barrier();
  auto end = current_time();
  transferTime += time_span(end, start);
  // Particles are now all stored on processor 0
}

void GFlowBase::recallParticlesByProcessor(vector<vector<Particle> >& allParticles) {
  // Get all the particles back from the sectorizations
  auto start = current_time();
  int max = min(ndx*ndy, numProc);
  if (rank==0) {
    vector<Particle> myParticles;
    for (auto &p : sectorization.getParticles())
      myParticles.push_back(p);
    allParticles.push_back(myParticles); // Zeroth entry
  }
  for (int proc=1; proc<max; ++proc) {
    if (rank==0) {
      int size = 0;
      // Recieve how many Particles to expect to recieve
      CommWork.Recv( &size, 1, MPI_INT, proc, 0);
      // Create a buffer of Particles to send
      Particle *buffer = new Particle[size];
      CommWork.Recv( buffer, size, PARTICLE, proc, 0);
      vector<Particle> myParticles;
      for (int i=0; i<size; ++i) myParticles.push_back(buffer[i]);
      allParticles.push_back(myParticles);
    }
    else if (rank==proc) {
      auto &parts = sectorization.getParticles();
      int size = parts.size(), root = 0;
      CommWork.Send( &size, 1, MPI_INT, root, 0);
      // Create a buffer of Particles to send
      Particle *buffer = new Particle[size];
      int i=0;
      for (auto p=parts.begin(); p!=parts.end(); ++p, ++i) buffer[i] = *p;
      CommWork.Send( buffer, size, PARTICLE, root, 0);
    }
  }
  CommWork.Barrier();
  auto end = current_time();
  transferTime += time_span(end, start);
  // Particles are now all stored, sorted by processor, on processor 0
}

double GFlowBase::reduceStatFunction(StatFunc f, int choice) {
  int total = 0;
  pair<double, int> data = sectorization.doStatFunction(f);
  if (choice==0) { // Average
    double send = data.first*data.second, aggregate = 0;
    CommWork.Reduce(&send, &aggregate, 1, MPI_DOUBLE, MPI_SUM, 0);
    CommWork.Reduce(&data.second, &total, 1, MPI_INT, MPI_SUM, 0);
    if (rank==0) return total>0 ? aggregate/total : 0;
    else return 0;
  }
  else if (choice==1) { // Max
    double aggregate;
    CommWork.Reduce(&data, &aggregate, 1, MPI_DOUBLE, MPI_MAX, 0);
    return aggregate;
  }
  else if (choice==2) { // SUM
    double aggregate;
    CommWork.Reduce(&data, &aggregate, 1, MPI_DOUBLE, MPI_SUM, 0);
    return aggregate;
  }
  else return -1;
}

list<Particle> GFlowBase::createParticles(vector<vec2> positions, double radius, double dispersion, std::function<vec2(double)> velocity, std::function<double(double)> omega, double coeff, double dissipation, double repulsion, int interaction) {
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (auto pos : positions) {
      double r = dispersion>0 ? (1-drand48()*dispersion)*radius : radius;
      double theta = drand48()*2*PI;
      Particle p(pos, r);
      p.theta = theta;
      p.velocity = velocity(p.invMass);
      p.omega = omega(p.invII);
      p.dissipation = dissipation;
      p.coeff = coeff;
      p.repulsion = repulsion;
      p.interaction = interaction;
      allParticles.push_back(p);
    }
  }
  return allParticles;
}

void GFlowBase::createAndDistributeParticles(int number, const Bounds &b, Sectorization &sectors, double radius, double dispersion, std::function<vec2(double)> velocity, double coeff, double dissipation, double repulsion, int interaction) {
  // Bounds width and height
  double width = b.right-b.left, height = b.top-b.bottom;
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (int i=0; i<number; ++i) {
      double r = dispersion>0 ? (1-drand48()*dispersion)*radius : radius;
      vec2 pos(drand48()*(width-2*r)+r, drand48()*(height-2*r)+r);
      Particle p(pos, r);
      p.velocity = velocity(p.invMass);
      p.dissipation = dissipation;
      p.coeff = coeff;
      p.repulsion = repulsion;
      allParticles.push_back(p);
    }
  }
  distributeParticles(allParticles, sectors);
}

void GFlowBase::createAndDistributeParticles(const vector<double>& radii, const vector<int>& interaction, const Bounds &b, Sectorization &sectors, std::function<vec2(double)> velocity, double coeff, double dissipation, double repulsion) {
  int number = min(interaction.size(), radii.size());
  // Bounds width and height
  double width = b.right-b.left, height = b.top-b.bottom;
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (int i=0; i<number; ++i) {
      double r = radii.at(i);
      vec2 pos(drand48()*(width-2*r)+r, drand48()*(height-2*r)+r);
      Particle p(pos, r);
      p.velocity = velocity(p.invMass);
      p.dissipation = dissipation;
      p.coeff = coeff;
      p.repulsion = repulsion;
      p.interaction = interaction.at(i);
      allParticles.push_back(p);
    }
  }
  distributeParticles(allParticles, sectors);
}

void GFlowBase::createAndDistributeParticles(const vector<vec2>& positions, const vector<double>& radii, const vector<int>& interaction, const Bounds &b, Sectorization &sectors, std::function<vec2(double)> velocity, double coeff, double dissipation, double repulsion) {
  int number = min(positions.size(), radii.size());
  // Bounds width and height
  double width = b.right-b.left, height = b.top-b.bottom;
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (int i=0; i<number; ++i) {
      double r = radii.at(i);
      Particle p(positions.at(i), r);
      p.velocity = velocity(p.invMass);
      p.dissipation = dissipation;
      p.coeff = coeff;
      p.repulsion = repulsion;
      p.interaction = interaction.at(i);
      allParticles.push_back(p);
    }
  }
  distributeParticles(allParticles, sectors);
}

// Data will be sent to the rank 0 (head) processor
vector<vec2> GFlowBase::findPackedSolution(int number, double radius, Bounds b, vec2 force, double expandTime, double relaxTime) {
  // Create a sectorization
  Sectorization packedSectors;
  packedSectors.setSimBounds(b);
  // Vectors for the four corners
  vec2 ll(b.left, b.bottom), lr(b.right, b.bottom), tl(b.left, b.top), tr(b.right, b.top);
  // Add the boundary walls
  packedSectors.addWall(Wall(ll,lr));
  packedSectors.addWall(Wall(lr,tr));
  packedSectors.addWall(Wall(tl,tr));
  packedSectors.addWall(Wall(ll,tl));
  // Add any "real" walls (so we don't put particles inside walls)
  for (auto &w : walls) packedSectors.addWall(w);
  // Calculate parameters
  double initialRadius = 0.2*radius, finalRadius = 1.*radius;
  int expandSteps = expandTime/epsilon, relaxSteps = relaxTime/epsilon;
  double dr = (finalRadius-initialRadius)/static_cast<double>(expandSteps);
  // Set up sectorization
  packedSectors.setASize(number); //** AD HOC
  packedSectors.setDoInteractions(doInteractions);
  setUpSectorization(packedSectors, force);
  packedSectors.setDrag(default_packed_drag);
  createAndDistributeParticles(number, b, packedSectors, radius, 0, ZeroV); //** START WITH INITIAL RADIUS, NOT RADIUS
  packedSectors.initialize();
  // Simulate motion and expansion
  int printIter = 0;
  double counter=0, delay=1./15.;
  for (int i=0; i<expandSteps; ++i) {
    for (auto &p : packedSectors.getParticles()) p.sigma += dr; //** DOESN'T DO ANYTHING
    packedSectors.update();
    // For observation
    if (writeCreation) {
      counter+=epsilon;
      if (delay<counter) {
	printToCSV(writeDirectory+"Init/fps"+toStr(printIter)+".csv", packedSectors.getParticles());
	++printIter;
	counter=0;
      }
    }
  }
  // Simulate pure motion (relaxation)
  for (int i=0; i<relaxSteps; ++i) {
    packedSectors.update();
    // For observation
    if (writeCreation) {
      counter+=epsilon;
      if (delay<counter) {
        printToCSV(writeDirectory+"Init/fps"+toStr(printIter)+".csv", packedSectors.getParticles());
        ++printIter;
        counter=0;
      }
    }
  }
  // Get positions
  auto sectorParticles = packedSectors.getParticles();
  if (numProc>1) {
    int size = sectorParticles.size();
    // Get the number of particles to expect from each processor
    int *sizeBuff = 0;
    if (rank==0) sizeBuff = new int[numProc];
    CommWork.Gather(&size, 1, MPI_INT, sizeBuff, 1, MPI_INT, 0);
    // Find the max number of particles in any sector with the head processor and broadcast it back
    int max=0;
    if (rank==0)
      for (int i=0; i<numProc; ++i)
	if (max<sizeBuff[i]) max = sizeBuff[i];
    CommWork.Bcast(&max, 1, MPI_INT, 0);
    // Now everyone knows what the max is. Allocate arrays only as large as neccessary
    Particle *parts = new Particle[max], *buffer = 0;
    if (rank==0) buffer = new Particle[number*max];
    int i=0;
    for (auto p=sectorParticles.begin(); p!=sectorParticles.end(); ++p, ++i) 
      parts[i] = *p;
    // Send Particle data to master processor
    CommWork.Gather(parts, size, PARTICLE, buffer, max, PARTICLE, 0);
    // If the head processor, fill with particle positions
    vector<vec2> positions;
    if (rank==0) {
      positions.reserve(number);
      for (int r=0; r<numProc; ++r)
	for (int i=0; i<sizeBuff[r]; ++i) {
	  int add = max*r+i;
	  positions.push_back(buffer[add].position);
	}
    }
    // Delete memory
    if (parts) delete[] parts;
    if (buffer) delete[] buffer;
    if (sizeBuff) delete[] sizeBuff;
  // Return positions
    return positions;
  }
  else { // Single processor run
    vector<vec2> positions;
    for (auto p : sectorParticles) positions.push_back(p.position);
    // Return positions
    return positions;
  }
}

vector<vec2> GFlowBase::findPackedSolution(const vector<double>& radii, const vector<int>& interactions, const Bounds& b, vec2 force, double expandTime, double relaxTime) {
  int number = min(radii.size(), interactions.size());
  // Create a sectorization
  Sectorization packedSectors;
  packedSectors.setSimBounds(b);
  // Vectors for the four corners
  vec2 ll(b.left, b.bottom), lr(b.right, b.bottom), tl(b.left, b.top), tr(b.right, b.top);
  // Add the boundary walls
  packedSectors.addWall(Wall(ll,lr));
  packedSectors.addWall(Wall(lr,tr));
  packedSectors.addWall(Wall(tl,tr));
  packedSectors.addWall(Wall(ll,tl));
  // Add any "real" walls (so we don't put particles inside walls)
  for (auto &w : walls) packedSectors.addWall(w);
  // Calculate parameters
  double initialFraction = 0.2, finalFraction = 1.;
  int expandSteps = expandTime/epsilon, relaxSteps = relaxTime/epsilon;
  double dr = (finalFraction-initialFraction)/static_cast<double>(expandSteps);
  double fraction = initialFraction;
  // Set up sectorization
  packedSectors.setASize(number); //** AD HOC
  packedSectors.setDoInteractions(doInteractions);
  setUpSectorization(packedSectors, force);
  packedSectors.setDrag(default_packed_drag);
  // Create random initial positions
  createAndDistributeParticles(radii, interactions, b, packedSectors);
  packedSectors.initialize();
  // Simulate motion and expansion
  for (int i=0; i<expandSteps; ++i) {
    int j=0;
    for (auto &p : packedSectors.getParticles()) {
      p.sigma = radii.at(j)*fraction;
      ++j;
    }
    fraction += dr;
    packedSectors.update();
  }
  // Simulate pure motion (relaxation)
  for (int i=0; i<relaxSteps; ++i) packedSectors.update();
  // Get positions
  auto sectorParticles = packedSectors.getParticles();
  if (numProc>1) {
    int size = sectorParticles.size();
    // Get the number of particles to expect from each processor
    int *sizeBuff = 0;
    if (rank==0) sizeBuff = new int[numProc];
    CommWork.Gather(&size, 1, MPI_INT, sizeBuff, 1, MPI_INT, 0);
    // Find the max number of particles in any sector with the head processor and broadcast it back
    int max=0;
    if (rank==0)
      for (int i=0; i<numProc; ++i)
	if (max<sizeBuff[i]) max = sizeBuff[i];
    CommWork.Bcast(&max, 1, MPI_INT, 0);
    // Now everyone knows what the max is. Allocate arrays only as large as neccessary
    Particle *parts = new Particle[max], *buffer = 0;
    if (rank==0) buffer = new Particle[number*max];
    int i=0;
    for (auto p=sectorParticles.begin(); p!=sectorParticles.end(); ++p, ++i)
      parts[i] = *p;
    // Send Particle data to master processor
    CommWork.Gather(parts, size, PARTICLE, buffer, max, PARTICLE, 0);
    // If the head processor, fill with particle positions
    vector<vec2> positions;
    if (rank==0) {
      positions.reserve(number);
      for (int r=0; r<numProc; ++r)
	for (int i=0; i<sizeBuff[r]; ++i) {
	  int add = max*r+i;
	  positions.push_back(buffer[add].position);
	}
    }
    // Delete memory
    if (parts) delete[] parts;
    if (buffer) delete[] buffer;
    if (sizeBuff) delete[] sizeBuff;
    // Return positions
    return positions;
  }
  else { // Single processor run
    vector<vec2> positions;
    for (auto p : sectorParticles) positions.push_back(p.position);
    // Return positions
    return positions;
  }
}

vector<vec2> GFlowBase::findLatticeSolution(int number, double radius, Bounds bounds, int ltype, double drf, double prob) {
  if (ltype<0) ltype = latticeType;
  // Calculate some parameters
  double latt = sqrt(3.);
  double dr = drf*radius;
  double width = bounds.right - bounds.left;
  int nx = width/(2.*radius+dr);
  int ny = (bounds.top-bounds.bottom)/(2.*radius+dr);
  
  int count = 0;
  vector<vec2> positions;
  if (ltype==0) { // Hexagonal lattice
    double dx = 0, Y = radius;
    for (int y=0; y<ny && count<number; ++y) {
      double X = radius+dr;
      if (y%2==1) X += (radius+0.5*dr); // Staggered
      for (int x=0; x<nx && X+radius<width; ++x) {
	if (prob<1. && drand48()<prob) {
          positions.push_back(vec2(X,Y));
          ++count;
        }
        else {
          positions.push_back(vec2(X,Y));
          ++count;
        }
	// Update
	X += (2.*radius+dr);
      }
      Y += (latt*radius+dr);
    }
  }
  else { // Rectangular lattice
    double Y = radius;
    for (int y=0; y<ny && count<number; ++y) {
      double X = radius+dr;
      for (int x=0; x<nx && count<number; ++x) {
	if (prob<1. && drand48()<prob) {
	  positions.push_back(vec2(X,Y));
	  ++count; 
	}
	else {
	  positions.push_back(vec2(X,Y));
	  ++count;
	}
	// Update x
	X += (2*radius+dr);
      }
      Y += (2*radius+dr);
    }
  }

  return positions;
}
  
// ----- TO GO TO GFLOW.CPP -----
  
void GFlowBase::setAsBacteria() {
  Characteristic *B = new Bacteria;
  sectorization.setCharacteristic(B);
  delete B;
  sectorization.setDrag(true);
  sectorization.stopParticles();
  sectorization.setASize(5000); //** AD HOC
  doFields = true;
  Bounds bounds(left, right, bottom, top);
  Resource.setBounds(bounds);
  Resource.setWrap(true);
  Resource.setResolution(0.025);
  Resource.setDiffusion(default_resource_diffusion);
  Resource.set(ConstantField); // Give some initial resource so we don't die immediately
  Resource.setLambda(default_resource_lambda);
  Waste.setBounds(bounds);
  Waste.setResolution(0.025);
  Waste.setWrap(true);
  Waste.setDiffusion(default_waste_diffusion);
  Waste.setLambda(default_waste_lambda);
  // Set termination condition
  addTerminationCondition(Stat_Number_Particles, allGone);
}

void GFlowBase::createSquare(int number, double radius, double width, double height, double vsgma, double dispersion, int interaction) {
  // Start the clock
  auto start = current_time();
  // Discard any old state
  discard();
  // Bounds
  Bounds bounds(0, width, 0, height);
  setBounds(bounds);
  // Everyone knows where the walls are
  Wall w[4];
  w[0] = Wall(left, bottom, right, bottom);
  w[1] = Wall(left, bottom, left, top);
  w[2] = Wall(left, top, right, top);
  w[3] = Wall(right, bottom, right, top);
  for (int i=0; i<4; ++i) {
    w[i].dissipation = 0;
    w[i].coeff = 0;
    addWall(w[i]);
  }
  gravity = 0; 
  // Velocity initialization function
  std::function<vec2(double)> velocity = [&] (double invMass) { 
    double angle = 2*PI*drand48();
    vec2 v(cos(angle), sin(angle));
    double ke = fabs(vsgma*randNormal());
    double velocity = sqrt(2*invMass*ke/127.324);
    v *= velocity;
    return v;
  };
  // Create 
  // Create particles and distribute them to the processors
  vector<vec2> positions = findPackedSolution(number, radius, bounds);  
  // Create particles at the given positions with - Radius, Dispersion, Velocity function, Angular velocity function, Coeff, Dissipation, Repulsion, Interaction
  particles = createParticles(positions, radius, dispersion, velocity, ZeroOm, 0, 0, default_sphere_repulsion, interaction);
  // Send out particles
  setUpSectorization();
  distributeParticles(particles, sectorization);
  sectorization.initialize();
  // sectorization.setASize(number); //**  AD HOC
  // End setup timing
  auto end = current_time();
  bubbleBounds = visBounds = NullBounds;
  setUpTime = time_span(end, start);
}

void GFlowBase::createBuoyancyBox(double radius, double bR, double density, double width, double depth, double velocity, double dispersion, int interaction) {
  // Start the clock
  auto start = current_time();
  bR = fabs(bR); radius = fabs(radius);
  // Discard any old state
  discard();
  // Bounds
  double height = depth+2*bR+10*radius;
  Bounds bounds(0, width, 0, height);
  Bounds initialBounds(0, width, 0, depth);
  setBounds(bounds);
  // Everyone knows where the walls are
  addWall(left, bottom, right, bottom);
  addWall(left, bottom, left, top);
  addWall(left, top, right, top);
  addWall(right, bottom, right, top);
  // Set up sectorization
  gravity = vec2(0.,-1.);
  setUpSectorization();
  // Calculate how many particles we should start with
  double maxPack = PI/(2*sqrt(3)); // Hexagonal packing
  double Vfill = 0.9*width*height;
  // Get positions
  vector<vec2> positions;
  vector<double> radii;
  vector<int> interactions;
  int number = 0;
  if (dispersion==0) { // All particles have the same radius
    double Vgrain = PI*sqr(radius);
    number = 0.85*maxPack*Vfill/Vgrain;
    // Create particles and distribute them to the processors
    if (latticeType>-1) positions = findLatticeSolution(number, radius, initialBounds, latticeType);
    else positions = findPackedSolution(number, radius, initialBounds, gravity);
  }
  else {
    // Not all particles have the same radius
    double totalV = 0;
    // Get enough particles to fill the volume
    while (totalV<Vfill) {
      double r = (1-drand48()*dispersion)*radius;
      totalV += PI*sqr(r);
      radii.push_back(r);
      interactions.push_back(interaction);
    }
    number = radii.size();
    if (latticeType>-1) positions = findLatticeSolution(number, radius, initialBounds, latticeType);
    else positions = findPackedSolution(radii, interactions, bounds, gravity);
  }
  // Create particles at the given positions with - Radius, Dispersion, Velocity function, Angular velocity function, Coeff, Dissipation, Repulsion, Interaction
  double coeff = default_sphere_coeff;
  if (dispersion==0) {
    particles = createParticles(positions, radius, dispersion, ZeroV, ZeroOm, coeff, default_sphere_dissipation, default_sphere_repulsion, interaction);
    distributeParticles(particles, sectorization);
  }
  else createAndDistributeParticles(positions, radii, interactions, bounds, sectorization, ZeroV, coeff);
  // Add the large ball
  sectorization.initialize();
  // sectorization.setASize(number+1); //** AD HOC
  // End setup timing
  auto end = current_time();
  bubbleBounds = visBounds = NullBounds;
  setUpTime = time_span(end, start);
}

bool GFlowBase::loadBuoyancy(string fileName, double radius, double velocity, double density, bool drag) {
  // Start the clock
  auto start = current_time();
  // Discard any old state
  discard();
  // Load data
  if (!loadConfigurationFromFile(fileName)) return false;
  // Lambda that finds the tops of the balls
  StatFunc upperEdge = [&] (const vector<Particle> &particles, int&) {
    double tp = -1e9;
    for (const auto &p : particles)
      if (tp<p.position.y+p.sigma) tp = p.position.y+p.sigma;
    return tp;
  };
  // Find the tops of the balls
  double tp = reduceStatFunction(upperEdge, 1);
  // Set bounds
  double height = tp + 10*radius+0.5;
  Bounds bounds(left, right, bottom, height);
  setBounds(bounds);
  // Set use characteristics
  sectorization.setUseCharacteristics(drag);
  // Clear old walls, create new ones
  walls.clear();
  addWall(left, bottom, right, bottom);
  addWall(left,bottom,left,top);
  addWall(right,bottom,right,top);
  addWall(left,top,right,top);
  // Send out particles
  setUpSectorization();
  sectorization.initialize();
  // Add the intruding particle
  if (0<radius) {
    Particle P((right-left)/2, tp+radius, radius);
    vec2 V(0, -velocity);
    P.velocity = V;
    P.setDensity(density);
    if (drag) sectorization.insertParticle(P, new ConstantVelocity(V));
    else sectorization.insertParticle(P);
  }
  // Set a termination condition
  addTerminationCondition(Stat_Large_Object_Height, belowZero);
  // End setup timing
  auto end = current_time();
  setUpTime = time_span(end, start);
  bubbleBounds = visBounds = NullBounds;
  return true;
}

void GFlowBase::addStatFunction(StatFunc sf, string str) {
  statFunctions.push_back(pair<StatFunc,string>(sf, str));
  statRecord.push_back(vector<vec2>());
}

void GFlowBase::addStatPlot(StatPlot sp, string str, double lower, double upper) {
  statPlots.push_back(pair<StatPlot,string>(sp, str));
  statPlotRecord.push_back(vector<double>(statPlotBins));
  statPlotBounds.push_back(pair<double,double>(lower,upper));
}

void GFlowBase::addTerminationCondition(StatFunc sf, std::function<bool(double)> tc) {
  terminationConditions.push_back(pair<StatFunc, std::function<bool(double)> >(sf, tc));
}

string GFlowBase::printStatFunctions(string label) {
  if (statFunctions.empty() || rank!=0) return "";
  stringstream stream;
  string str, strh;
  for (int i=0; i<statFunctions.size(); ++i) {
    string name = statFunctions.at(i).second;
    stream << name << label << "=" << mmPreproc(statRecord.at(i));
    stream >> strh;
    str += (strh+";\nPrint[\""+name+"\"]\nListLinePlot["+name+label+",PlotStyle->Black,ImageSize->Large,PlotRange->All]\n");
    stream.clear(); strh.clear();
  }
  return str;
}

string GFlowBase::printStatPlots(string label) {
  if (statPlots.empty() || rank!=0) return "";
  stringstream stream;
  string str, strh;
  for (int i=0; i<statPlots.size(); ++i) {
    string name = statPlots.at(i).second;
    stream << name << label << "=" << mmPreproc(statPlotRecord.at(i));
    stream >> strh;
    str += (strh+";\nPrint[\""+name+"\"]\nListLinePlot["+name+label+",PlotStyle->Black,ImageSize->Large,PlotRange->All]\n");
    stream.clear(); strh.clear();
  }
  return str;
}

vector<vpair> GFlowBase::getWallsPositions() {
  vector<vpair> positions;
  for (auto w : walls) positions.push_back(vpair(w.getLeft(), w.getRight()));
  return positions;
}

string GFlowBase::printWallsCommand() {
  stringstream stream;
  string str, strh, range;
  // Create range
  stream << "{{" << left << "," << right << "},{" << bottom << "," << top << "}}";
  stream >> range;
  stream.clear();
  // Create walls graphic
  stream << "w=" << mmPreproc(getWallsPositions(),3) << ";";
  stream >> str;
  stream.clear();
  str += "\n";
  if (!walls.empty()) {
    stream << "walls=Graphics[{";
    for (int i=0; i<walls.size(); ++i) {
      stream << "{Blue,Thick,Line[{w[[" << i+1 << "]][[1]],w[[" << i+1 << "]][[2]]}]}";
      if (i!=walls.size()-1) stream << ",";
    }
    stream << "},PlotRange->" + range + "];\n";
  }
  else stream << "walls={};\n";
  stream >> strh;
  return str+strh;
}

string GFlowBase::printAnimationCommand(int mode, bool novid, string label) {
  stringstream stream;
  string command, strh, range, scl;
  
  stream << "pos" << label << "=" << printPositionRecord(mode) << ";";
  stream >> command;
  stream.clear();
  command += "\n";

  // Print animation bounds if neccessary
  if (!animationBounds.empty()) {
    stream << "bnds=" << mmPreproc(animationBounds,2) << ";";
    stream >> strh;
    stream.clear();
    command += (strh+"\n");
  }

  stream << "len=" << recIter << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "scale=" << scale << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  // Triangle animation command
  command += "tri[dt_]:=Triangle[{dt[[1]]+dt[[2]]*{Cos[dt[[3]]],Sin[dt[[3]]]},dt[[1]]+dt[[2]]*{Cos[dt[[3]]+2*Pi/3],Sin[dt[[3]]+2*Pi/3]},dt[[1]]+dt[[2]]*{Cos[dt[[3]] + 4*Pi/3], Sin[dt[[3]] + 4*Pi/3]}}];\n";
  // Disk animation command
  command += "dsk[tr_]:={Black,Disk[tr[[1]],tr[[2]]]};\n";
  // Oriented Disk animation command
  command += "odsk[tr_]:={{Black,{Disk[tr[[1]],tr[[2]]]}},{Red,Line[{tr[[1]], tr[[1]] + tr[[2]]{Cos[tr[[3]]],Sin[tr[[3]]]}}]}};\n";
  // Point animation command
  command += "pnt[tr_]:={Black,Point[tr[[1]]]};\n";
  // Dot (point) animation command
  command += "dot[tr_]:={Black,Point[tr]};\n";

  // Bounds boxes
  int animationCentering=0; //**
  if (animationCentering==1) {
    command += "Bd[{a_,b_,c_,d_}]:={{a,c},{a,d},{b,d},{b,c}};\n";
    command += "bounds=Table[rect[Bd[bnds[[i]]]],{i,1,Length[bnds]}];\n";
    command += "Bd2[{a_,b_,c_,d_}]:={{a,b},{c,d}};\n";
  }

  // Create range
  if (animationCentering==1) stream << "Bd2[bnds[[i]]]";
  else stream << "{{" << left << ","<< right << "},{" << bottom << "," << top << "}}";
  stream >> range;
  stream.clear();

  stream << "ImageSize->{";
  if (animationCentering==1) stream << "scale*(bnds[[i]][[2]]-bnds[[i]][[1]]),scale*(bnds[[i]][[4]]-bnds[[i]][[3]])";
  else stream << "scale*" << right-left << ",scale*" << top-bottom;
  stream << "}";
  stream >> scl;
  stream.clear();
  // Print walls
  command += printWallsCommand();

  stream << (mode==0 ? "dsk" : "dot");
  stream >> strh;

  command += "disks=Table[Graphics[Table[" + strh + "[pos" + label + "[[i]][[j]]],{j,1,Length[pos" + label + "[[i]]]}],PlotRange->" + range + "],{i,1,len}];\n";
  command += ("frames=Table[Show[disks[[i]],walls," + scl + "],{i,1,len}];\n");
  if (!novid) command += "Export[\"vid.avi\",frames,\"CompressionLevel\"->0];\n";
  command += "ListAnimate[frames]";
  
  return command;  
}

string GFlowBase::printSpecialAnimationCommand(bool novid) {
  stringstream stream;
  string command, strh, range, scl;

  stream << "pos=" << mmPreproc(specialRecord,3) << ";";
  stream >> command;
  stream.clear();
  command += "\n";

  stream << "len=" << recIter << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "scale=" <<scale << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "{{" << left << "," << right << "},{" << bottom << "," << top << "}}";
  stream >> range;
  stream.clear();

  stream << "ImageSize->{scale*" << right-left << ",scale*" << top-bottom << "}";
  stream >> scl;
  stream.clear();
  // Print walls
  command += printWallsCommand();
  command += "nproc=";
  stream << numProc;
  stream >> strh;
  stream.clear();
  command += strh+";\n";
  command += "colors = Table[RGBColor[RandomReal[],RandomReal[],RandomReal[]],{i,1,nproc}]\n";
  command += "sdisk[tr_]:={colors[[tr[[3]]+1]],Disk[tr[[1]],tr[[2]]]};\n";
  command += "disks=Table[ Graphics[ Table[sdisk[pos[[i]][[j]] ],{j,1,Length[pos[[i]]]}],PlotRange->" + range + "], {i,1,len}];\n";
  command += ("frames=Table[Show[disks[[i]],walls," + scl + "],{i,1,len}];\n");
  if (!novid) command += "Export[\"vid.avi\",frames,\"CompressionLevel\"->0];\n";
  command += "ListAnimate[frames]";

  return command;
}

string GFlowBase::printForcesAnimationCommand(int mode, bool novid) {
  stringstream stream;
  string command, strh, range, scl;
  // Find maximum force
  double maxF = 0;
  for (const auto& v : positionRecord)
    for (auto f : v) {
      double force = std::get<4>(f);
      if (force>maxF) maxF = force;
    }
  // Create command - do this as the mode switching for now
  if (mode==0) stream << "pos=" << mmPreproc(positionRecord,3) << ";"; 
  else stream << "pos=" << mmPreproc(positionRecord) << ";";
  stream >> command;
  stream.clear();
  command += "\n";
  // Print max force
  stream << "maxF=" << maxF << ";";
  stream >> strh;
  stream.clear();
  command += (strh+'\n');
  // Print length and scale
  stream << "len=" << recIter << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "scale=" <<scale << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  // Create range string
  stream << "{{" << left << "," << right << "},{" << bottom << "," << top << "}}";
  stream >> range;
  stream.clear();
  // Create scale string
  stream << "ImageSize->{scale*" << right-left << ",scale*" << top-bottom << "}";
  stream >> scl;
  stream.clear();
  // Create walls graphic
  command += printWallsCommand();
  // Create disk graphics
  command += "col[f_]:=RGBColor[f/maxF,0,0];\n";
  command += "sdisk[tr_]:={col[tr[[3]]],Disk[tr[[1]],tr[[2]]]};\n";
  command += "disks=Table[Graphics[Table[sdisk[pos[[i]][[j]] ],{j,1,Length[pos[[i]]]}],PlotRange->" + range + "], {i,1,len}];\n";
  command += ("frames=Table[Show[disks[[i]],walls," + scl + "],{i,1,len}];\n");
  if (!novid) command += "Export[\"vid.avi\",frames,\"CompressionLevel\"->0];\n";
  command += "ListAnimate[frames]";

  return command;
}

string GFlowBase::printBulkAnimationCommand(bool novid, bool center) {
  stringstream stream;
  string command, strh;
  command = "scale=" + toStr(scale) + ";\n";
  // Print bulk data
  stream << "bulk=" << mmPreproc(bulkRecord,2) << ";";
  stream >> strh;
  stream.clear();
  command += (strh + "\n");
  // Print bounds data
  stream << "bnds=" << mmPreproc(bulkBounds,2) << ";";
  stream >> strh;
  stream.clear();
  command += (strh + "\n");
  // Bounds boxes
  command += "rect[b_]:=Graphics[{Opacity[0],EdgeForm[Green],Polygon[b]}];\n";
  command += "Bd[{a_,b_,c_,d_}]:={{a,c},{a,d},{b,d},{b,c}};\n";
  command += "bounds=Table[rect[Bd[bnds[[i]]]],{i,1,Length[bnds]}];\n";
  if (center) command += "Bd2[{a_,b_,c_,d_}]:={{a,b},{c,d}};\n";
  // Add walls
  if (!center) command += printWallsCommand();
  strh.clear();
  stream << "bulkFrames=Table[Show[{";
  if (!center) stream << "walls,bounds[[i]],";
  stream << "Graphics[{Thick,Line[bulk[[i]]]}]},PlotRange->";
  if (center) stream << "Bd2[bnds[[i]]]";
  else stream << "{{" << left << ","<< right << "},{" << bottom << "," << top << "}}";
  stream << ",ImageSize->{";
  if (center) stream << "scale*(bnds[[i]][[2]]-bnds[[i]][[1]]),scale*(bnds[[i]][[4]]-bnds[[i]][[3]])";
  else stream << "scale*" << right-left << ",scale*" << top-bottom;
  stream << "}],{i,1,Length[bulk]}];";
  stream >> strh;
  if (!novid) strh += "Export[\"vidB.avi\",bulkFrames,\"CompressionLevel\"->0];\n";
  strh += "ListAnimate[bulkFrames]";
  return command+strh;
}

string GFlowBase::printSnapshot() {
  stringstream stream;
  string command, strh, range, scl;

  vector<Particle> allParticles;
  recallParticles(allParticles);
  int number = allParticles.size();
  // Find average radius
  double radius = 0.;
  for (auto &p : allParticles) radius += p.sigma;
  radius = number>0 ? radius/number : 0.;
  // Record positions
  vector<vec2> positions;
  for (const auto &p : allParticles) positions.push_back(p.position);
  stream << "snap=" << mmPreproc(positions,2) << ";\n";
  stream >> command;
  stream.clear();
  command += "\n";
  stream << "R=" << radius << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "len=" << positions.size() << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "scale=" <<scale << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  // Triangle animation command
  //command += "tri[dt_]:=Triangle[{dt[[1]]+dt[[2]]*{Cos[dt[[3]]],Sin[dt[[3]]]},dt[[1]]+dt[[2]]*{Cos[dt[[3]]+2*Pi/3],Sin[dt[[3]]+2*Pi/3]},dt[[1]]+dt[[2]]*{Cos[dt[[3]] + 4*Pi/3],Sin[dt[[3]] + 4*Pi/3]}}];\n";
  // Disk animation command
  command += "dsk[tr_]:={Black,Disk[tr,R]};\n";
  // Oriented Disk animation command
  //command += "odsk[tr_]:={{Black,{Disk[tr[[1]],tr[[2]]]}},{Red,Line[{tr[[1]],tr[[1]]+tr[[2]]{Cos[tr[[3]]],Sin[tr[[3]]]}}]}};\n";
  // Point animation command
  //command += "pnt[tr_]:={Black,Point[tr[[1]]]};\n";
  // Dot (point) animation command
  //command += "dot[tr_]:={Black,Point[tr]};\n";
  // Create range
  stream << "{{" << left << "," << right << "},{" << bottom << "," << top << "}}";
  stream >> range;
  stream.clear();

  stream << "ImageSize->{scale*" << right-left << ",scale*" << top-bottom << "}";
  stream >> scl;
  stream.clear();

  // Print walls
  command += printWallsCommand();
  command += "disks=Graphics[Table[dsk[snap[[i]]], {i,1,len} ],PlotRange->" + range + "];\n";
  command += ("img=Show[disks,walls," + scl + "]\n");
  return command;
}

void GFlowBase::printSectors() {
  string command;
  for (int y=0; y<ndy; ++y)
    for (int x=0; x<ndx; ++x) {
      if (rank == ndx*y+x) {
        cout << "sec" << x << "x" << y << "=" << sectorization.printSectors() << ";\n";
        cout << "MatrixPlot[sec" << x << "x" << y << "]\n";
      }
      CommWork.Barrier();
    }
}

void GFlowBase::getBulkData(vector<VPair>& lines, double volCutoff, double boxV, double dr, double upperVolCutoff) {
  vector<Particle> allParticles;
  recallParticles(allParticles);
  getBulkData(allParticles, lines, Bounds(left, right, bottom, top), volCutoff, boxV, dr, upperVolCutoff);
}

vector<double> GFlowBase::getBulkData(vector<Particle> &allParticles, Bounds region, double volCutoff, double boxV, double dr, double upperVolCutoff) {
  string str="-1";
  vector<VPair> lines;
  return getBulkData(allParticles, str, lines, false, region, volCutoff, boxV, dr, upperVolCutoff);
}

vector<double> GFlowBase::getBulkData(vector<Particle> &allParticles, vector<VPair>& lines, Bounds region, double volCutoff, double boxV, double dr, double upperVolCutoff) {
  string str="-1";
  return getBulkData(allParticles, str, lines, true, region, volCutoff, boxV, dr, upperVolCutoff);
}

vector<double> GFlowBase::getBulkData(vector<Particle> &allParticles, string &shapes, vector<VPair>& lines, bool getOutline, Bounds region, double volCutoff, double boxV, double dr, double upperVolCutoff) {
  // Might use full bounds
  if (region.right<region.left || region.top<region.bottom) region = Bounds(left, right, bottom, top);
  // Invalid bounds
  if (region.right<=region.left || region.top<=region.bottom) return vector<double>();
  // Calculate parameters
  double sx = sqrt(boxV), sy = sx;
  int nsx = (region.right-region.left)/sx, nsy = (region.top-region.bottom)/sy;
  sx = (region.right-region.left)/nsx; sy = (region.top-region.bottom)/nsy;
  ++nsx; ++nsy;
  list<int> *sectors = new list<int>[nsx*nsy];
  // Fill sectors
  int i=0;
  for (const auto &p : allParticles) {
    int sec_x = (p.position.x - region.left)/sx;
    int sec_y = (p.position.y - region.bottom)/sy;
    if (-1<sec_x && sec_x<nsx && -1<sec_y && sec_y<nsy)
      sectors[nsx*sec_y+sec_x].push_back(i);
    // If part of the particle sticks into the region, place the particle in the closest sector. This ignores particles "diagonal" to the region that stick into the region
    else if (region.left<p.position.x+p.sigma   && -1<sec_y && sec_y<nsy)
      sectors[nsx*sec_y+0].push_back(i); // Left
    else if (p.position.x-p.sigma<region.right  && -1<sec_y && sec_y<nsy)
      sectors[nsx*sec_y+(nsx-1)].push_back(i); // Right
    else if (region.bottom<p.position.y+p.sigma && -1<sec_x && sec_x<nsx)
      sectors[sec_x].push_back(i); // Bottom
    else if (p.position.y-p.sigma<region.top    && -1<sec_x && sec_x<nsx)
      sectors[nsx*(nsy-1)+sec_x].push_back(i); // Top
    ++i;
  }
  // Find the particle of maximum radius
  double maxR = 0;
  for (const auto &p : allParticles)
    if (p.sigma>maxR) maxR = p.sigma;
  // Check which sector centers are covered by particles
  int *array = new int[nsx*nsy];
  int sweepX = (maxR+dr)/sx, sweepY = (maxR+dr)/sy;
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x) {
      // Check whether center of sector is covered by any particles
      array[nsx*y+x] = nsx*y+x;
      bool done = false;
      vec2 pos((x+0.5)*sx+region.left, (y+0.5)*sy+region.bottom);
      // If outside the simulation bounds, there are no bubbles here
      if (!Bounds(left,right,bottom,top).contains(pos)) {
	array[nsx*y+x] = -1;
	continue;
      }
      // Sweep through sectors
      int startX = max(0, x-sweepX), endX = min(nsx-1, x+sweepX);
      int startY = max(0, y-sweepY), endY = min(nsy-1, y+sweepY);
      // Check your own sector first for speed's sake
      for (const auto index : sectors[nsx*y+x]) {
	vec2 dR = pos-allParticles.at(index).position;
	double R = dr+allParticles.at(index).sigma;
	if (sqr(dR)<sqr(R)) {
	  done = true;
	  array[nsx*y+x] = -1;
	  break;
	}
      }
      // Search the other sectors
      for (int Y=startY; Y<endY && !done; ++Y) {
	for (int X=startX; X<endX && !done; ++X) {
	  if (X!=x || Y!=y)
	    for (const auto index : sectors[nsx*Y+X]) {
	      vec2 dR= pos-allParticles.at(index).position;
	      double R = dr+allParticles.at(index).sigma;
	      if (sqr(dR)<sqr(R)) {		
		done = true;
		array[nsx*y+x] = -1;
		break;
	      }
	    }
	  if (done) break;
	}
	if (done) break;
      }
    }
  // Unite bubbles
  unite(array, nsx, nsy);
  // Find all bubbles that border an edge
  std::set<int> edgeBubbles;
  auto contains = [&] (int i) { return edgeBubbles.find(i)!=edgeBubbles.end(); };
  auto add = [&] (int i) {
    int h = getHead(array, i);
    if (h<0 || contains(h));
    else edgeBubbles.insert(h);
  };
  for (int x=0; x<nsx; ++x) add(x); // Bottom
  for (int x=0; x<nsx; ++x) add(nsx*nsy-x-1); // Top
  for (int y=0; y<nsy; ++y) add(nsx*y); // Left
  for (int y=0; y<nsy; ++y) add(nsx*y+nsx-1); // Right
  // Point all sectors to their head
  for (int k=0; k<nsx*nsy; ++k) array[k] = getHead(array,k);
  // "Erase" bubbles that touch the bounds
  for (int k=0; k<nsx*nsy; ++k)
    if (contains(array[k])) array[k] = -1;
  // Collect all head nodes
  std::map<int, int> labels; // Maps old label to new label
  std::map<int, int> volCount; // Maps (new) label to number of cells in the bubble
  int lab = 0;
  for (int k=0; k<nsx*nsy; ++k)
    if (array[k]!=-1)
      // If we have not already recorded this head
      if (labels.find(array[k])==labels.end()) {
	labels.insert(pair<int, int>(array[k], lab));
	volCount.insert(pair<int, int>(lab, 0));
	++lab;
      }
  // If there are no empty volumes
  if (labels.empty()) {
      delete [] sectors;
      delete [] array;
      return vector<double>();
  }
  // Count volumes
  for (int k=0; k<nsx*nsy; ++k) {
    if (array[k]!=-1) { // Search for the proper index
      int j=0;
      // Increment volume counter
      auto it = volCount.find(array[k]);
      if (it!=volCount.end())
	++it->second;
    }
  }
  // Record volume sizes
  vector<double> bubbles;
  i=0;
  for (const auto c : volCount) {
    double vol = sx*sy*c.second;
    if (volCutoff<vol && vol<upperVolCutoff) bubbles.push_back(vol);
  }
  // Remove volumes that are to small
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x) {
      int j = array[nsx*y+x];
      if (-1<j) {
	auto it = volCount.find(j);
	if (sx*sy*it->second<volCutoff) array[nsx*y+x] = -1;
      }
    }
  // Record data in scalar field
  if (!bubbles.empty() && bubbleField.empty()) {
    // Set up field
    if (bubbleField.empty()) {
      bubbleField.setBounds(region);
      bubbleField.setResolution(sx);
      bubbleField.setPrintPoints(200);
    }
    // Update
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x)
	bubbleField.at(x,y) += (array[nsx*y+x]>-1 ? 1 : 0);
  }
  // Create outline
  if (getOutline) createOutline(array, nsx, nsy, sx, sy, region, lines);
  // Record the picture of the volumes as a matrix
  if (shapes!="-1") createMatrix(array, nsx, nsy, sx, sy, volCutoff, volCount, shapes);

  // Sort sizes
  std::sort(bubbles.begin(), bubbles.end());
  // Clean up and return
  delete [] sectors;
  delete [] array;
  return bubbles;
}

inline void GFlowBase::unite(int *array, int nsx, int nsy) {
  // Find out which volumes are connected
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x) {
      // Check the sectors around you, point yourself to the largest head
      int h0=-1, h1=-1, h2=-1, h3=-1;
      if (array[nsx*y+x]!=-1) {
        int head = getHead(array, nsx*y+x);
        if (0<x && -1<array[nsx*y+x-1]) { // Left
	  h0 = getHead(array, nsx*y+(x-1));
	  if (h0<head) head = h0;
	}
        if (x+1<nsx && -1<array[nsx*y+x+1]) { // Right
	  h1 = getHead(array, nsx*y+(x+1));
	  if (h1<head) head = h1;
	}
        if (0<y && -1<array[nsx*(y-1)+x]) { // Bottom
	  h2 = getHead(array, nsx*(y-1)+x);
	  if (h2<head) head = h2;
	}
        if (y+1<nsy && -1<array[nsx*(y+1)+x]) { // Top
	  h3 = getHead(array, nsx*(y+1)+x);
	  if (h3<head) head = h3;
	}
	// Set your head and the heads of all your neighbors as well
	array[nsx*y+x] = head;
	if (-1<h0) array[h0] = head;
	if (-1<h1) array[h1] = head;
	if (-1<h2) array[h2] = head;
	if (-1<h3) array[h3] = head;
      }
    }
}

inline int GFlowBase::getHead(int* array, int index) {
  if (index<0) return -1;
  while (array[index]!=index && -1<array[index]) index = array[index];
  return index;
}

inline void GFlowBase::createOutline(int *array, int nsx, int nsy, double sx, double sy, Bounds region, vector<VPair>& lines) {
  // Create a lambda function to edge detect
  auto isEdge = [&] (int x, int y) {
    if (x<0 || nsx<=x || y<0 || nsy<=y || array[nsx*y+x]!=-1) return false;
    if (x+1<nsx && array[nsx*y+x+1]>-1)  return true;
    else if (0<x && array[nsx*y+x-1]>-1) return true;
    else if (y+1<nsy && array[nsx*(y+1)+x]>-1) return true;
    else if (0<y && array[nsx*(y-1)+x]>-1) return true;
    return false;
  };
  lines.clear();
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x) {
      if (isEdge(x,y)) {
	vec2 pos((x+0.5)*sx+region.left, (y+0.5)*sy+region.bottom);
	if (isEdge(x+1,y-1)) lines.push_back(VPair(pos, pos+vec2(sx,-sy)));
	if (isEdge(x+1,y))   lines.push_back(VPair(pos, pos+vec2(sx,0)));
	if (isEdge(x+1,y+1)) lines.push_back(VPair(pos, pos-vec2(sx,sy)));
	if (isEdge(x,y+1))   lines.push_back(VPair(pos, pos+vec2(0,sy)));
      }
    }
}

inline void GFlowBase::createMatrix(int* array, int nsx, int nsy, double sx, double sy, double volCutoff, std::map<int,int> volCount, string& shapes) {
  stringstream stream;
  stream << "{";
  for (int y=nsy-1; 0<=y; --y) {
    stream << "{";
    for (int x=0; x<nsx; ++x) {
      if (-1<array[nsx*y+x] && (volCount.empty() || volCutoff<volCount.find(array[nsx*y+x])->second*sx*sy))
	stream << array[nsx*y+x]+1;
      else stream << -1;
      if (x!=nsx-1) stream << ",";
    }
    stream << "}";
    if (y!=0) stream << ",";
  }
  stream << "}";
  stream >> shapes;
}

string GFlowBase::printPositionRecord(int mode) {
  switch (mode) {
  default:
  case 0: { // Normal
    return mmPreproc(positionRecord, 3);
    break;
  }
  case 1: { // Compressed form
    vector<vector<vec2> > reducedData;
    for (const auto& v : positionRecord) {
      vector<vec2> positions;
      for (const auto& p : v) positions.push_back(std::get<0>(p));
      reducedData.push_back(positions);
    }
    return mmPreproc(reducedData, 2);
  }
  case 2: { // Long form
    return mmPreproc(positionRecord);
  }
  }
}

string GFlowBase::printResource() {
  stringstream stream;
  string str;
  stream << Resource;
  stream >> str;
  return str;
}

string GFlowBase::printWaste() {
  stringstream stream;
  string str;
  stream << Waste;
  stream >> str;
  return str;
}

inline Bounds GFlowBase::followBallBounds() {
  Bounds domain;
  double X = reduceStatFunction(Stat_Large_Object_X);
  double Y = reduceStatFunction(Stat_Large_Object_Height);
  static double R = reduceStatFunction(Stat_Large_Object_Radius, 1);
  domain.bottom = Y-2*R;         domain.top = Y+18*R;
  domain.left = max(left,X-8*R); domain.right = min(right, X+8*R);
  return domain;
}

inline void GFlowBase::updateFitness() {
  if (Fitness.empty()) {
    Fitness.setBounds(Resource.getBounds());
    Fitness.setResolution(Resource.getResolution());
  }
  int nsx = Fitness.getNSX(), nsy = Fitness.getNSY();
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x) {
      double res = Resource.at(x,y), wst = Waste.at(x,y);
      double fitness = alphaR*res/(res+csatR) - alphaW*wst/(wst+csatW);
      Fitness.at(x,y) = fitness;
    }
}
