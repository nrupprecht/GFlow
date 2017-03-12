#include "GFlowBase.h"

using MPI::COMM_WORLD;

GFlowBase::GFlowBase() {
  left = 0; right = 1; bottom = 0; top = 1;
  wrapX = true; wrapY = true;
  gravity = vect<>(0., -1.);
  temperature = 0; viscosity = 1.308e-3;
  time = 0;
  epsilon = 1e-4; sqrtEpsilon = sqrt(epsilon);
  sectorization.setEpsilon(epsilon);
  dispTime = 1./15.; lastDisp = 0;
  iter = 0; recIter = 0; maxIter = 0;
  runTime = 0;
  running = false;
  transferTime = 0;

  doInteractions = true;

  doWork = false;
  cutoff = 0.05;
  skinDepth = 0.25*cutoff;

  //---
  setUpTime = 0;
  recPositions = false;
  recKE = false;
  //---

  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();

  // Define the particle data type
  MPI_Type_contiguous( 16, MPI_DOUBLE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );
}

GFlowBase::~GFlowBase() {}

void GFlowBase::initialize() {
  // Rank 0 distributes, other processes listen and obey
  
}

void GFlowBase::run(double runLength) {
  resetVariables();
  // setUpSectorization();
  // Calculate the number of iterations we will run for
  maxIter = runLength/epsilon;
  clock_t start = clock();
  running = true;
  record(); // Initial record
  for (iter=0; iter<maxIter; ++iter) {
    if (doWork) {
      objectUpdates();
    }
    if (time-lastDisp>dispTime) record();

    logisticUpdates();
    COMM_WORLD.Barrier();
  }
  running = false;
  clock_t end = clock();
  runTime = (double)(end-start)/CLOCKS_PER_SEC;
  COMM_WORLD.Barrier();
  gatherData();
  COMM_WORLD.Barrier();
}

void GFlowBase::addWall(Wall w) {
  walls.push_back(w);
  sectorization.addWall(w);
}

void GFlowBase::addWall(double sx, double sy, double ex, double ey) {
  Wall w(vect<>(sx,sy), vect<>(ex,ey));
  walls.push_back(w);
  sectorization.addWall(w);
}

void GFlowBase::addParticle(Particle p) {
  particles.push_back(p);
}

void GFlowBase::addParticle(double x, double y, double r) {
  particles.push_back(Particle(x, y, r));
}

void GFlowBase::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
  sectorization.setSimBounds(l, r, b, t);
}

void GFlowBase::setBounds(Bounds b) {
  left = b.left; right = b.right; bottom = b.bottom; top = b.top;
  sectorization.setSimBounds(b);
}

bool GFlowBase::loadConfigurationFromFile (string filename) {
  throw false; // UNIMPLEMENTED
}

bool GFlowBase::createConfigurationFile (string filename) {
  throw false; // UNIMPLEMENTED
}

void GFlowBase::setUpSectorization() {
  // Decide how to divide up the space into domains
  bestProcessorGrid(ndx, ndy, numProc, Bounds(left, right, bottom, top));
  // Calculate what bounds this processor is in charge of
  if (rank<ndx*ndy) {
    doWork = true; // This processor needs to do work
    sectorization.setSimBounds(left, right, bottom, top);
    Bounds domainBounds = getBoundsForProc(rank);
    sectorization.setBounds(domainBounds);
    sectorization.setGravity(gravity);
    sectorization.setCutoff(cutoff);
    sectorization.setSkinDepth(skinDepth);
    sectorization.setDoInteractions(doInteractions);
    sectorization.initialize();
  }
  else doWork = false;
}

void GFlowBase::setUpSectorization(Sectorization &sectors, double cutoff, double skinDepth) {
  // Decide how to divide up the space into domains
  bestProcessorGrid(ndx, ndy, numProc, sectors.getSimBounds() );
  // Calculate what bounds this processor is in charge of
  if (rank<ndx*ndy) {
    doWork = true; // This processor needs to do work
    sectors.setSimBounds(left, right, bottom, top);
    Bounds domainBounds = getBoundsForProc(rank, sectors.getSimBounds() );
    sectors.setBounds(domainBounds);
    sectors.setGravity(gravity);
    sectors.setCutoff(cutoff);
    sectors.setSkinDepth(skinDepth);
    sectors.setDoInteractions(doInteractions);
    sectors.initialize();
  }
  else doWork = false;
}

void GFlowBase::resetVariables () {
  time = 0;
  iter = 0;
  recIter = 0;
  maxIter = 0;
  runTime = 0;
}

void GFlowBase::objectUpdates () {
  sectorization.update();
}

void GFlowBase::logisticUpdates() {
  time += epsilon;
}

void GFlowBase::record() {
  if (recPositions || recKE) recordPositions();
  // Update display 
  lastDisp = time;
  ++recIter;
}

void GFlowBase::resets() {

}

void GFlowBase::gatherData() {
  
}

void GFlowBase::discard() {
  walls.clear();
  particles.clear();
  sectorization.discard();
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
    for (auto p=allParticles.begin(); p!=allParticles.end(); ++p) {
      if (bnds.contains(p->position)) {
        remove.push_back(p);
        parts.push_back(*p);
      }
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
      COMM_WORLD.Send( &size, 1, MPI_INT, proc, 0);
      // Send the actual data to processor i
      COMM_WORLD.Send( buffer, size, PARTICLE, proc, 0);
    }
    else if (rank==proc) { // Recieve
      int size = 0, root = 0;
      // Recieve the amount of data we should expect
      COMM_WORLD.Recv( &size, 1, MPI_INT, root, 0);
      // Recieve the actual data
      Particle *buffer = new Particle[size];
      COMM_WORLD.Recv( buffer, size, PARTICLE, root, 0);
      // Add particles to sectors
      for (int i=0; i<size; ++i) sectors.addParticle(buffer[i]);
    }
  }
  COMM_WORLD.Barrier();
}

list<Particle> GFlowBase::createParticles(vector<vect<> > positions, double radius, double dispersion, std::function<vect<>(double)> velocity, double coeff, double dissipation, double repulsion, int interaction) {
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (auto pos : positions) {
      double r = dispersion>0 ? (1-drand48()*dispersion)*radius : radius;
      Particle p(pos, r);
      p.velocity = velocity(p.invMass);
      p.dissipation = dissipation;
      p.coeff = coeff;
      p.repulsion = repulsion;
      allParticles.push_back(p);
    }
  }
  return allParticles;
}

void GFlowBase::createAndDistributeParticles(int number, const Bounds &b, Sectorization &sectors, double radius, double dispersion, std::function<vect<>(double)> velocity, double coeff, double dissipation, double repulsion, int interaction) {
  // Bounds width and height
  double width = b.right-b.left, height = b.top-b.bottom;
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (int i=0; i<number; ++i) {
      double r = dispersion>0 ? (1-drand48()*dispersion)*radius : radius;
      vect<> pos(drand48()*(width-2*r)+r, drand48()*(height-2*r)+r);
      Particle p(pos, r);
      p.velocity = velocity(p.invMass);
      p.dissipation = dissipation;
      p.coeff = coeff;
      p.repulsion = repulsion;
      allParticles.push_back(p);
    }
  }

  distributeParticles(allParticles, sectors);

  /*
  // Take your own particles
  if (rank==0) {
    Bounds bnds = getBoundsForProc(0, sectors.getSimBounds() );
    vector<list<Particle>::iterator> remove;
    for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
      if (bnds.contains(p->position)) {
        remove.push_back(p);
        particles.push_back(*p);
      }
    for (auto &p : remove) allParticles.erase(p);
    // Add particles to sectors
    for (auto &p : particles) sectors.addParticle(p);
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
      COMM_WORLD.Send( &size, 1, MPI_INT, proc, 0);
      // Send the actual data to processor i
      COMM_WORLD.Send( buffer, size, PARTICLE, proc, 0);
    }
    else if (rank==proc) { // Recieve
      int size = 0, root = 0;
      // Recieve the amount of data we should expect
      COMM_WORLD.Recv( &size, 1, MPI_INT, root, 0);
      // Recieve the actual data
      Particle *buffer = new Particle[size];
      COMM_WORLD.Recv( buffer, size, PARTICLE, root, 0);
      // Add particles to sectors
      for (int i=0; i<size; ++i) sectors.addParticle(buffer[i]);
    }
  }
  COMM_WORLD.Barrier();
  */
}

// Data will be sent to the rank 0 (head) processor
vector<vect<> > GFlowBase::findPackedSolution(int number, double radius, Bounds b, vect<> force, double expandTime, double relaxTime) {
  // Create a sectorization
  Sectorization packedSectors;
  packedSectors.setSimBounds(b);
  // Vectors for the four corners
  vect<> ll(b.left, b.bottom), lr(b.right, b.bottom), tl(b.left, b.top), tr(b.right, b.top);
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
  // Distribute particles
  createAndDistributeParticles(number, b, packedSectors, radius);
  // Set up sectorization
  setUpSectorization(packedSectors, 2*radius, 0.25*radius);
  packedSectors.setGravity(force); // Have to do this after setUpSectorization
  packedSectors.setDrag(default_packed_drag);
  // Simulate motion and expansion
  for (int i=0; i<expandSteps; ++i) {
    for (auto &p : packedSectors.getParticles()) p.sigma += dr;
    packedSectors.update();
  }
  // Simulate pure motion (relaxation)
  for (int i=0; i<relaxSteps; ++i) packedSectors.update();
  // Get positions
  auto sectorParticles = packedSectors.getParticles();
  int size = sectorParticles.size();
  Particle *parts = new Particle[size], *buffer = 0;
  int i=0;
  for (auto p=sectorParticles.begin(); p!=sectorParticles.end(); ++p, ++i) parts[i] = *p;
  // If we didn't know that the total number of particles is [number], we would have to use MPI reduce operation
  if (rank==0) buffer = new Particle[number];
  // Send data to master processor
  COMM_WORLD.Gather(parts, size, PARTICLE, buffer, number, PARTICLE, 0);
  if (rank==0) {
    vector<vect<> > positions;
    positions.reserve(number);
    for (int i=0; i<number; ++i) positions.push_back(buffer[i].position);
    return positions;
  }
  else return vector<vect<> >(); // Empty vector
}

// ----- TO GO TO GFLOW.CPP -----

void GFlowBase::createSquare(int number, double radius, double width, double height, double vsgma, double dispersion) {
  // Start the clock
  auto begin = clock();
  // Discard any old state
  discard();
  // Bounds
  Bounds bounds(0, width, 0, height);
  setBounds(bounds);
  // Set cutoff and skin depth
  cutoff = 2*radius;
  skinDepth = 0.25*radius;
  // Everyone knows where the walls are
  addWall(left, bottom, right, bottom);
  addWall(left, bottom, left, top);
  addWall(left, top, right, top);
  addWall(right, bottom, right, top);
  // Set up sectorization
  setUpSectorization();
  gravity = 0;
  sectorization.setGravity(gravity);
  // Velocity initialization function
  std::function<vect<>(double)> velocity = [&] (double invMass) { 
    double angle = 2*PI*drand48();
    vect<> v(cos(angle), sin(angle));
    double ke = fabs(vsgma*randNormal());
    double velocity = sqrt(2*invMass*ke/127.324);
    v *= velocity;
    return v;
  };
  // Create particles and distribute them to the processors
  vector<vect<> > positions = findPackedSolution(number, radius, bounds);  
  list<Particle> allParticles = createParticles(positions, radius, dispersion, velocity, 0, 0);
  // Send out particles
  distributeParticles(allParticles, sectorization);
  // End setup timing
  auto end = clock();
  setUpTime = (double)(end-begin)/CLOCKS_PER_SEC;
}

void GFlowBase::createBuoyancyBox(double radius, double bR, double density, double width, double depth, double velocity, double dispersion) {
  // Start the clock
  auto begin = clock();
  // Discard any old state
  discard();
  // Bounds
  Bounds bounds(0, width, 0, depth+2*bR);
  setBounds(bounds);
  // Set cutoff and skin depth
  cutoff = radius+bR;
  skinDepth = 0.25*radius;
  // Everyone knows where the walls are
  addWall(left, bottom, right, bottom);
  addWall(left, bottom, left, top);
  addWall(left, top, right, top);
  addWall(right, bottom, right, top);
  // Set up sectorization
  setUpSectorization();
  gravity = vect<>(0,-1);
  sectorization.setGravity(gravity);
  // Calculate how many particles we should start with
  double maxPack = PI/(2*sqrt(3)); // Hexagonal packing
  double Vfill = width*depth, Vgrain = PI*sqr(radius*(1-0.5*dispersion));
  int number = maxPack * Vfill / Vgrain;
  // Create particles and distribute them to the processors
  vector<vect<> > positions = findPackedSolution(number, radius, bounds, gravity);
  // Only allow particles that are below depth
  list<list<Particle>::iterator> remove;
  list<Particle> allParticles = createParticles(positions, radius, dispersion, ZeroV);
  for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
    if (depth<p->position.y+p->sigma) remove.push_back(p);
  for (auto &p : remove) allParticles.erase(p);
  // Add the large ball
  if (rank==0) {
    Particle ball(width/2, depth+bR, bR);
    ball.setDensity(density);
    allParticles.push_back(ball);
  }
  // Send out particles
  distributeParticles(allParticles, sectorization);
  // End setup timing
  auto end = clock();
  setUpTime = (double)(end-begin)/CLOCKS_PER_SEC;
}

void GFlowBase::recordPositions() {
  // Get all the particles back from the sectorizations
  auto begin = clock();
  vector<Particle> allParticles;
  int max = min(ndx*ndy, numProc);
  if (rank==0)
    for (auto &p : sectorization.getParticles())
      allParticles.push_back(p);
  for (int proc=1; proc<max; ++proc) {
    if (rank==0) {
      int size = 0;
      // Recieve how many Particles to expect to recieve
      COMM_WORLD.Recv( &size, 1, MPI_INT, proc, 0);
      // Create a buffer of Particles to send
      Particle *buffer = new Particle[size];
      COMM_WORLD.Recv( buffer, size, PARTICLE, proc, 0);
      for (int i=0; i<size; ++i) allParticles.push_back(buffer[i]);
    }
    else if (rank==proc) {
      list<Particle> &parts = sectorization.getParticles();
      int size = parts.size(), root = 0;
      COMM_WORLD.Send( &size, 1, MPI_INT, root, 0);
      // Create a buffer of Particles to send
      Particle *buffer = new Particle[size];
      int i=0;
      for (auto p=parts.begin(); p!=parts.end(); ++p, ++i) buffer[i] = *p;
      COMM_WORLD.Send( buffer, size, PARTICLE, root, 0);
    }
  }
  COMM_WORLD.Barrier();
  auto end = clock();
  transferTime += (double)(end-begin)/CLOCKS_PER_SEC;
  // Particles are now all stored on processor 0
  
  if (recPositions) {
    vector<pair<vect<>, double> > positions;
    for (auto p : allParticles) positions.push_back(pair<vect<>, double>(p.position, p.sigma));
    positionRecord.push_back(positions);
  }
  if (recKE) {
    double ke = 0;
    for (auto p : allParticles) ke +=sqr(p.velocity)/p.invMass;
    ke *= (0.5*(1./allParticles.size()));
    keRecord.push_back(ke);
  }
}

vector<vpair> GFlowBase::getWallsPositions() {
  vector<vpair> positions;
  for (auto w : walls) positions.push_back(vpair(w.left, w.right));
  return positions;
}

string GFlowBase::printAnimationCommand(bool novid) {
  stringstream stream;
  string command, strh, range, scale;
  
  stream << "pos=" << mmPreproc(positionRecord) << ";";
  stream >> command;
  stream.clear();
  command += "\n";
  stream << "w=" << mmPreproc(getWallsPositions()) << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "len=" << recIter << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\nscale=100;\n");

  stream << "{{" << left << "," << right << "},{" << bottom << "," << top << "}}";
  stream >> range;
  stream.clear();

  stream << "ImageSize->{scale*" << right-left << ",scale*" << top-bottom << "}";
  stream >> scale;
  stream.clear();

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
  stream.clear();
  command += strh;
  command += "disk[tr_]:={Black,Disk[tr[[1]],tr[[2]]]};\n";
  command += "disks=Table[ Graphics[ Table[disk[pos[[i]][[j]] ],{j,1,Length[pos[[i]]]}],PlotRange->" + range + "], {i,1,len}];\n";
  command += ("frames=Table[Show[disks[[i]],walls," + scale + "],{i,1,len}];\n");
  if (!novid) command += "Export[\"vid.avi\",frames,\"CompressionLevel\"->0];\n";
  command += "ListAnimate[frames]";
  
  return command;  
}
