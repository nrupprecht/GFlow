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
  dispTime = 1./15.; lastDisp = -2*dispTime;
  iter = 0; recIter = 0; maxIter = 0;
  runTime = 0;
  running = false;

  doInteractions = true;

  //---
  recPositions = false;
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
  setUpSectorization();
  // Calculate the number of iterations we will run for
  maxIter = runLength/epsilon;
  clock_t start = clock();
  running = true;

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

void GFlowBase::addWall (Wall w) {
  walls.push_back(w);
}

void GFlowBase::addParticle (Particle p) {
  particles.push_back(p);
}

void GFlowBase::addParticle (double x, double y, double r) {
  particles.push_back(Particle(x, y, r));
}

void GFlowBase::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
  needsRemake = true;
  sectorization.setSimBounds(l, r, b, t);
}

void GFlowBase::setBounds(Bounds b) {
  left = b.left; right = b.right; bottom = b.bottom; top = b.top;
  needsRemake = true;
  sectorization.setSimBounds(b);
}

bool GFlowBase::loadConfigurationFromFile (string filename) {
  throw false; // UNIMPLEMENTED
}

bool GFlowBase::createConfigurationFile (string filename) {
  throw false; // UNIMPLEMENTED
}

void GFlowBase::setUpSectorization () {
  // Decide how to divide up the space into domains
  //-- FIGURE OUT HOW TO DO THIS BEST LATER
  ndx = 4; ndy = 4;
  
  // Calculate what bounds this processor is in charge of
  if (rank<ndx*ndy) {
    doWork = true; // This processor needs to do work
    sectorization.setSimBounds(left, right, bottom, top);
    Bounds domainBounds = getBoundsForProc(rank);
    sectorization.setBounds(domainBounds);
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
  if (recPositions) recordPositions();
  // Update display 
  lastDisp = time;
  ++recIter;
}

void GFlowBase::resets() {

}

void GFlowBase::gatherData() {
  
}

Bounds GFlowBase::getBoundsForProc(int rnk) {
  double dx = (right-left)/ndx, dy = (top-bottom)/ndy;
  double l = (rnk%ndx)*dx, r = l+dx;
  double b = (rnk/ndx)*dy, t = b+dy;
  return Bounds(l, r, b, t);
}

// ----- TO GO TO GFLOW.CPP -----

void GFlowBase::createSquare(int N, double radius, double width, double height, double vsgma) {
  Bounds bounds(0, width, 0, height);
  setBounds(bounds);
  setUpSectorization();
  // Set cutoff and skin depth
  cutoff = 2*radius;
  skinDepth = 0.25*radius;
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (int i=0; i<N; ++i) {
      vect<> pos(drand48()*width, drand48()*height);
      double angle = 2*PI*drand48();
      vect<> v(cos(angle), sin(angle));
      v *= (vsgma*randNormal());
      Particle p(pos, radius);
      p.velocity = v;
      allParticles.push_back(p);
    }
  }
  // Take your own particles
  if (rank==0) {
    Bounds b = getBoundsForProc(0);
    vector<list<Particle>::iterator> remove;
    for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
      if (b.contains(p->position)) {
        remove.push_back(p);
        particles.push_back(*p);
      }
    for (auto &p : remove) allParticles.erase(p);
    // Add particles to sectors
    for (auto &p : particles) sectorization.addParticle(p);
  }
  // Distribute particles to processes
  int max = min(ndx*ndy, numProc);
  for (int proc=1; proc<max; ++proc) {
    if (rank==0) { // Send
      vector<Particle> domainParticles;
      Bounds b = getBoundsForProc(proc);
      vector<list<Particle>::iterator> remove;
      for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
        if (b.contains(p->position)) {
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
      for (int i=0; i<size; ++i) sectorization.addParticle(buffer[i]);
    }
  }
  COMM_WORLD.Barrier();
}

void GFlowBase::recordPositions() {
  // Get all the particles back from the sectorizations
  vector<Particle> allParticles;
  int max = min(ndx*ndy, numProc);
  if (rank==0)
    for (auto &p : sectorization.getParticles())
      allParticles.push_back(p);
  for (int proc=1; proc<max; ++proc) {
    if (rank==0) {
      int size = 0;
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
      for (auto p=parts.begin(); p!=parts.end(); ++p, ++i) {
        buffer[i] = *p;
      }
      COMM_WORLD.Send( buffer, size, PARTICLE, root, 0);
    }
  }
  COMM_WORLD.Barrier();
  // Particles are now all stored on processor 0
  
  vector<pair<vect<>, double> > positions;
  for (auto p : allParticles) positions.push_back(pair<vect<>, double>(p.position, p.sigma));
  positionRecord.push_back(positions);
}

string GFlowBase::printAnimationCommand() {
  stringstream stream;
  string command, strh;
  
  stream << "len=" << recIter << ";";
  stream >> command;
  stream.clear();

  //stream << "pos=" << positionRecord << ";";
  //stream >> strh;
  //stream.clear();
  command += (strh+"\nscale=100;\n");
  
  command += "disk[tr_]:={Black,Disk[tr[[1]],tr[[2]] ]};\n";
  command += "frames=Table[ Graphics[ Table[disk[pos[[i]][[j]] ],{j,1,Length[pos[[i]] ]}],PlotRange->";

  stream << "{{" << left << "," << right << "},{" << bottom << "," << top << "}}";
  stream >> strh;
  command += (strh+" ], {i,1,len}];\nExport[\"vid.avi\",frames,\"CompressionLevel\"->0];ListAnimate[frames]");

  return command;
   
}
