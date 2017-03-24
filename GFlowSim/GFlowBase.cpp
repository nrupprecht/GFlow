#include "GFlowBase.h"

using MPI::COMM_WORLD;

GFlowBase::GFlowBase() {
  left = 0; right = 1; bottom = 0; top = 1;
  wrapX = false; wrapY = false;
  gravity = vec2(0., -1.);
  temperature = 0; viscosity = 1.308e-3;
  time = 0;
  epsilon = 1e-4; sqrtEpsilon = sqrt(epsilon);
  sectorization.setEpsilon(epsilon);
  dispTime = 1./15.; lastDisp = 0;
  startRec = 0;
  iter = 0; recIter = 0; maxIter = 0;
  runTime = 0;
  running = false;
  transferTime = 0;

  doInteractions = true;

  doWork = false;
  cutoff = 0.05;
  skinDepth = cutoff;

  //---
  setUpTime = 0;
  recPositions = false;
  recKE = false;
  //---

  // Get MPI system data
  rank = COMM_WORLD.Get_rank();
  numProc = COMM_WORLD.Get_size();
  CommWork = COMM_WORLD;

  // Define the particle data type
  MPI_Type_contiguous( 16*sizeof(floatType)/sizeof(float), MPI_FLOAT, &PARTICLE );
  MPI_Type_commit( &PARTICLE );
}

GFlowBase::~GFlowBase() {}

void GFlowBase::run(double runLength) {
  resetVariables();
  // Create work communicator
  //** int color = doWork ? 1 : 0;
  //** MPI::Intercomm CommWork = COMM_WORLD.Split(color, rank);
  sectorization.setCommWork(CommWork);
  // Calculate the number of iterations we will run for
  maxIter = runLength/epsilon;
  clock_t start = clock();
  running = true;
  if (startRec<=0) record(); // Initial record
  if (doWork) {
    for (iter=0; iter<maxIter; ++iter) {
      objectUpdates();
      if (time-lastDisp>dispTime && startRec<time) record();
      
      logisticUpdates();
      CommWork.Barrier();
    }
  }
  running = false;
  clock_t end = clock();
  runTime = (double)(end-start)/CLOCKS_PER_SEC;
  // Release work comm and set sectorization comm back to comm world
  sectorization.resetComm();
  // CommWork.Free(); 
  // Update transfer time
  transferTime += sectorization.getTransferTime();
}

void GFlowBase::addWall(Wall w) {
  walls.push_back(w);
  sectorization.addWall(w);
}

void GFlowBase::addWall(floatType sx, floatType sy, floatType ex, floatType ey) {
  Wall w(vec2(sx,sy), vec2(ex,ey));
  walls.push_back(w);
  sectorization.addWall(w);
}

void GFlowBase::addParticle(Particle p) {
  particles.push_back(p);
}

void GFlowBase::addParticle(floatType x, floatType y, floatType r) {
  particles.push_back(Particle(x, y, r));
}

int GFlowBase::getSize() {
  int sz = sectorization.getSize();
  int size = 0;
  CommWork.Reduce(&sz, &size, 1, MPI_INT, MPI_SUM, 0);
  return size;
}

void GFlowBase::setBounds(floatType l, floatType r, floatType b, floatType t) {
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

void GFlowBase::setTemperature(floatType t) {
  temperature = t;
  sectorization.setTemperature(t);
}

void GFlowBase::setViscosity(floatType h) {
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

bool GFlowBase::loadConfigurationFromFile (string filename) {
  throw false; // UNIMPLEMENTED
}

bool GFlowBase::createConfigurationFile (string filename) {
  throw false; // UNIMPLEMENTED
}

void GFlowBase::setUpSectorization() {
  // Decide how to divide up the space into domains
  bestProcessorGrid(ndx, ndy, numProc, Bounds(left, right, bottom, top));
  sectorization.giveDomainInfo(ndx, ndy);
  // Calculate what bounds this processor is in charge of
  if (rank<ndx*ndy) {
    doWork = true; // This processor needs to do work
    sectorization.setSimBounds(left, right, bottom, top);
    Bounds domainBounds = getBoundsForProc(rank);
    sectorization.setASize(particles.size()); //**
    sectorization.setBounds(domainBounds);
    sectorization.setGravity(gravity);
    sectorization.setCutoff(cutoff);
    sectorization.setSkinDepth(skinDepth);
    sectorization.setDoInteractions(doInteractions);
    // sectorization.initialize(); // We will do this later
  }
  else doWork = false;
}

void GFlowBase::setUpSectorization(Sectorization &sectors, floatType cutoff, floatType skinDepth) {
  // Decide how to divide up the space into domains
  bestProcessorGrid(ndx, ndy, numProc, sectors.getSimBounds() );
  sectors.giveDomainInfo(ndx, ndy);
  // Calculate what bounds this processor is in charge of
  if (rank<ndx*ndy) {
    doWork = true; // This processor needs to do work
    sectors.setSimBounds(left, right, bottom, top);
    Bounds domainBounds = getBoundsForProc(rank, sectors.getSimBounds() );
    sectors.setASize(particles.size()); //**
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
  if (recPositions || recKE || !statFunctions.empty()) {
    // Get all the particles back from the sectorizations
    vector<Particle> allParticles;
    recallParticles(allParticles);
    // Record positions
    if (recPositions) {
      vector<PData> positions;
      for (const auto &p : allParticles) positions.push_back(PData(p.position, p.sigma, p.theta, p.interaction, 0));
      positionRecord.push_back(positions);
    }

    // Record stat function statistics
    int i=0;
    for (auto &sf : statFunctions) {
      statRecord.at(i).push_back(vec2(time, sf.first(allParticles)));
      ++i;
    }
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
  floatType rmin=1e9, r=0, l=(b.right-b.left)/(b.top-b.bottom);
  // Find the domain grid that is most square and uses the most number of processors
  for (int i=1; i<number; ++i) {
    int j=number/i;
    floatType factor = static_cast<floatType>(j)/static_cast<floatType>(i)*l;
    r = factor + 1./factor;
    if (r<rmin) { rmin = r; xm = i; ym = j; }
  } 
  x = xm; y = ym;
}

Bounds GFlowBase::getBoundsForProc(int rnk) {
  floatType dx = (right-left)/ndx, dy = (top-bottom)/ndy;
  floatType l = (rnk%ndx)*dx, r = l+dx;
  floatType b = (rnk/ndx)*dy, t = b+dy;
  return Bounds(l, r, b, t);
}

Bounds GFlowBase::getBoundsForProc(int rnk, const Bounds &bnds) {
  floatType dx = (bnds.right-bnds.left)/ndx, dy = (bnds.top-bnds.bottom)/ndy;
  floatType l = (rnk%ndx)*dx, r = l+dx;
  floatType b = (rnk/ndx)*dy, t = b+dy;
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
  auto begin = clock();
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
  auto end = clock();
  transferTime += (double)(end-begin)/CLOCKS_PER_SEC;
  // Particles are now all stored on processor 0
}

void GFlowBase::recallParticlesByProcessor(vector<vector<Particle> >& allParticles) {
  // Get all the particles back from the sectorizations
  auto begin = clock();
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
  auto end = clock();
  transferTime += (double)(end-begin)/CLOCKS_PER_SEC;
  // Particles are now all stored, sorted by processor, on processor 0
}

list<Particle> GFlowBase::createParticles(vector<vec2> positions, floatType radius, floatType dispersion, std::function<vec2(floatType)> velocity, floatType coeff, floatType dissipation, floatType repulsion, int interaction) {
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (auto pos : positions) {
      floatType r = dispersion>0 ? (1-drand48()*dispersion)*radius : radius;
      Particle p(pos, r);
      p.velocity = velocity(p.invMass);
      p.dissipation = dissipation;
      p.coeff = coeff;
      p.repulsion = repulsion;
      p.interaction = interaction;
      allParticles.push_back(p);
    }
  }
  return allParticles;
}

void GFlowBase::createAndDistributeParticles(int number, const Bounds &b, Sectorization &sectors, floatType radius, floatType dispersion, std::function<vec2(floatType)> velocity, floatType coeff, floatType dissipation, floatType repulsion, int interaction) {
  // Bounds width and height
  floatType width = b.right-b.left, height = b.top-b.bottom;
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (int i=0; i<number; ++i) {
      floatType r = dispersion>0 ? (1-drand48()*dispersion)*radius : radius;
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

// Data will be sent to the rank 0 (head) processor
vector<vec2> GFlowBase::findPackedSolution(int number, floatType radius, Bounds b, vec2 force, floatType expandTime, floatType relaxTime) {
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
  floatType initialRadius = 0.2*radius, finalRadius = 1.*radius;
  int expandSteps = expandTime/epsilon, relaxSteps = relaxTime/epsilon;
  floatType dr = (finalRadius-initialRadius)/static_cast<floatType>(expandSteps);
  // Distribute particles
  createAndDistributeParticles(number, b, packedSectors, radius, 0, ZeroV);
  // Set up sectorization
  setUpSectorization(packedSectors, 2*radius, 0.25*radius);
  packedSectors.setGravity(force); // Have to do this after setUpSectorization
  packedSectors.setDrag(default_packed_drag);
  packedSectors.setASize(number); //**
  packedSectors.initialize();
  // Simulate motion and expansion
  for (int i=0; i<expandSteps; ++i) {
    for (auto &p : packedSectors.getParticles()) p.sigma += dr;
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

// ----- TO GO TO GFLOW.CPP -----

void GFlowBase::createSquare(int number, floatType radius, floatType width, floatType height, floatType vsgma, floatType dispersion) {
  // Start the clock
  auto begin = clock();
  // Discard any old state
  discard();
  // Bounds
  Bounds bounds(0, width, 0, height);
  setBounds(bounds);
  // Set cutoff and skin depth
  cutoff = 2*radius;
  skinDepth = radius;
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
  
  // Set up sectorization
  setUpSectorization();
  gravity = 0; 
  sectorization.setGravity(gravity);
  // Velocity initialization function
  std::function<vec2(floatType)> velocity = [&] (floatType invMass) { 
    floatType angle = 2*PI*drand48();
    vec2 v(cos(angle), sin(angle));
    floatType ke = fabs(vsgma*randNormal());
    floatType velocity = sqrt(2*invMass*ke/127.324);
    v *= velocity;
    return v;
  };

  // Create particles and distribute them to the processors
  vector<vec2> positions = findPackedSolution(number, radius, bounds);  
  // Create particles at the given positions with - Radius, Dispersion, Velocity, Coeff, Dissipation, Repulsion, Interaction
  list<Particle> allParticles = createParticles(positions, radius, dispersion, velocity, default_sphere_coeff, 0);
  // Send out particles
  distributeParticles(allParticles, sectorization);
  sectorization.initialize();
  // End setup timing
  auto end = clock();
  setUpTime = (double)(end-begin)/CLOCKS_PER_SEC;
}

void GFlowBase::createBuoyancyBox(floatType radius, floatType bR, floatType density, floatType width, floatType depth, floatType velocity, floatType dispersion) {
  // Start the clock
  auto begin = clock();
  // Discard any old state
  discard();
  // Bounds
  Bounds bounds(0, width, 0, depth+2*bR+radius);
  setBounds(bounds);
  // Set cutoff and skin depth
  cutoff = radius+max(bR,radius);
  skinDepth = 0.5*radius; //**
  // Everyone knows where the walls are
  addWall(left, bottom, right, bottom);
  addWall(left, bottom, left, top);
  addWall(left, top, right, top);
  addWall(right, bottom, right, top);
  // Set up sectorization
  setUpSectorization();
  gravity = vec2(0,-1);
  sectorization.setGravity(gravity);
  // Calculate how many particles we should start with

  floatType maxPack = PI/(2*sqrt(3)); // Hexagonal packing
  floatType Vfill = width*depth, Vgrain = PI*sqr(radius*(1-0.5*dispersion));
  int number = 0.9 * maxPack * Vfill / Vgrain;
  sectorization.setASize(number); //** There should be a more efficient way to do this, memory-wise
  // Create particles and distribute them to the processors
  vector<vec2> positions = findPackedSolution(number, radius, bounds, 0);
  // Only allow particles that are below depth
  list<Particle> allParticles = createParticles(positions, radius, dispersion, ZeroV, 0);
  // Remove particles whose centers are above the line
  for (auto p=allParticles.begin(); p!=allParticles.end();) {
    if (depth<p->position.y+p->sigma) p = allParticles.erase(p);
    else ++p;
  }
  // Add the large ball
  if (rank==0 && bR>0) {
    Particle ball(width/2, depth+bR, bR);
    ball.setDensity(density);
    allParticles.push_back(ball);
  }
  // Send out particles
  distributeParticles(allParticles, sectorization);
  sectorization.initialize();
  // End setup timing
  auto end = clock();
  setUpTime = (double)(end-begin)/CLOCKS_PER_SEC;
}

void GFlowBase::addStatFunction(StatFunc sf, string str) {
  statFunctions.push_back(pair<StatFunc,string>(sf, str));
  statRecord.push_back(vector<vec2>());
}

string GFlowBase::printStatFunctions() {
  if (statFunctions.empty() || rank!=0) return "";
  stringstream stream;
  string str, strh;
  for (int i=0; i<statFunctions.size(); ++i) {
    string name = statFunctions.at(i).second;
    stream << name << "=" << mmPreproc(statRecord.at(i));
    stream >> strh;
    str += (strh+";\nPrint[\""+name+"\"]\nListLinePlot["+name+",PlotStyle->Black,ImageSize->Large,PlotRange->All]\n");
    stream.clear(); strh.clear();
  }
  return str;
}

vector<vpair> GFlowBase::getWallsPositions() {
  vector<vpair> positions;
  for (auto w : walls) positions.push_back(vpair(w.getLeft(), w.getRight()));
  return positions;
}

string GFlowBase::printAnimationCommand(bool novid) {
  stringstream stream;
  string command, strh, range, scale;
  
  stream << "pos=" << mmPreproc(positionRecord,3) << ";";
  stream >> command;
  stream.clear();
  command += "\n";
  stream << "w=" << mmPreproc(getWallsPositions(),3) << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\n");

  stream << "len=" << recIter << ";";
  stream >> strh;
  stream.clear();
  command += (strh+"\nscale=100;\n");

  // Triangle animation command
  command += "tri[dt_]:=Triangle[{dt[[1]]+dt[[2]]*{Cos[dt[[3]]],Sin[dt[[3]]]},dt[[1]]+dt[[2]]*{Cos[dt[[3]]+2*Pi/3],Sin[dt[[3]]+2*Pi/3]},dt[[1]]+dt[[2]]*{Cos[dt[[3]] + 4*Pi/3], Sin[dt[[3]] + 4*Pi/3]}}];\n";
  // Disk animation command
  command += "dsk[tr_]:={Black,Disk[tr[[1]],tr[[2]]]};\n";
  // Oriented Disk animation command
  command += "odsk[tr_]:={{Black,{Disk[tr[[1]],tr[[2]]]}},{Red,Line[{tr[[1]], tr[[1]] + tr[[2]]{Cos[tr[[3]]],Sin[tr[[3]]]}}]}};\n";
  // Create range
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

  command += "disks=Table[Graphics[Table[dsk[pos[[i]][[j]]],{j,1,Length[pos[[i]]]}],PlotRange->" + range + "],{i,1,len}];\n";
  command += ("frames=Table[Show[disks[[i]],walls," + scale + "],{i,1,len}];\n");
  if (!novid) command += "Export[\"vid.avi\",frames,\"CompressionLevel\"->0];\n";
  command += "ListAnimate[frames]";
  
  return command;  
}

string GFlowBase::printSpecialAnimationCommand(bool novid) {
  stringstream stream;
  string command, strh, range, scale;

  stream << "pos=" << mmPreproc(specialRecord) << ";";
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
  command += "nproc=";
  stream << numProc;
  stream >> strh;
  stream.clear();
  command += strh+";\n";
  command += "colors = Table[RGBColor[RandomReal[],RandomReal[],RandomReal[]],{i,1,nproc}]\n";
  command += "sdisk[tr_]:={colors[[tr[[3]]+1]],Disk[tr[[1]],tr[[2]]]};\n";
  command += "disks=Table[ Graphics[ Table[sdisk[pos[[i]][[j]] ],{j,1,Length[pos[[i]]]}],PlotRange->" + range + "], {i,1,len}];\n";
  command += ("frames=Table[Show[disks[[i]],walls," + scale + "],{i,1,len}];\n");
  if (!novid) command += "Export[\"vid.avi\",frames,\"CompressionLevel\"->0];\n";
  command += "ListAnimate[frames]";

  return command;
}
