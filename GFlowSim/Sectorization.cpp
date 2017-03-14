#include "Sectorization.h"

using MPI::COMM_WORLD;

Sectorization::Sectorization() : nsx(3), nsy(3), secWidth(0), secHeight(0), time(0), epsilon(1e-4), sqrtEpsilon(1e-4), transferTime(0), wrapX(false), wrapY(false), doInteractions(true), drag(0), gravity(Zero), temperature(0), viscosity(1.308e-3), tempDelay(5e-3), sqrtTempDelay(sqrt(5e-3)), lastTemp(0), particles(0), sectors(0), doWallNeighbors(true), cutoff(0.1), skinDepth(0.025), itersSinceBuild(0), buildDelay(20), numProc(1), rank(0) {
  // Get our rank and the number of processors
  rank =    COMM_WORLD.Get_rank();
  numProc = COMM_WORLD.Get_size();  
  // Define the particle data type
  int size = sizeof(Particle)/sizeof(MPI_DOUBLE); // Should be 16
  MPI_Type_contiguous( size, MPI_DOUBLE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );
}

Sectorization::~Sectorization() {
  if (sectors) delete [] sectors;
  if (particles) delete [] particles;
  sectors = 0;
}

void Sectorization::initialize() {
  // First estimate
  secWidth = secHeight = cutoff+skinDepth;
  nsx = (bounds.right-bounds.left)/secWidth; 
  nsy = (bounds.top-bounds.bottom)/secHeight;
  // Actual width and height
  secWidth = (bounds.right-bounds.left)/nsx; 
  secHeight = (bounds.top-bounds.bottom)/nsy; 
  // Add for edge sectors
  nsx += 2; nsy += 2;
  // Remake sectors
  if (sectors) delete [] sectors;
  sectors = new list<Particle*>[nsx*nsy+1]; // +1 For special sector
  // Remake particles
  if (particles) delete [] particles;
  size = plist.size();
  particles = new Particle[size];
  int i=0;
  for (auto &p : plist) {
    particles[i] = p;
    ++i;
  }
  for (int i=0; i<size; ++i) add(&particles[i]);

  // Create the initial neighborlist
  createNeighborLists();
  if (doWallNeighbors) createWallNeighborList();
  itersSinceBuild = 0;

  // Reset times
  time = 0;
  transferTime = 0;
  lastTemp = 0;
}

list<Particle>& Sectorization::getParticles() {
  updatePList();
  return plist;
}

void Sectorization::setBounds(double l, double r, double b, double t) {
  bounds = Bounds(l, r, b, t);
}

void Sectorization::setBounds(Bounds b) {
  bounds = b;
}

void Sectorization::setSimBounds(double l, double r, double b, double t) {
  simBounds = Bounds(l, r, b, t);
}

void Sectorization::setInteractionType(int inter) {
  for (int i=0; i<size; ++i) particles[i].interaction = inter;
}

void Sectorization::discard() {
  if (particles) delete [] particles;
  particles = 0;
  size = 0;
  if (sectors) for (int i=0; i<nsx*nsy+1; ++i) sectors[i].clear();
  neighborList.clear();
  wallNeighbors.clear();
  walls.clear();
}

void Sectorization::setSimBounds(Bounds b) {
  simBounds = b;
}

void Sectorization::particleInteractions() {
  if (!doInteractions) return;
  // Neighbor list
  double n, s;
  for (auto &nl : neighborList) {
    auto p = nl.begin(); // The particle whose list this is
    auto q = p; ++q;
    for (; q!=nl.end(); ++q) { // Try to interact with all other particles in the nl
      vect<> displacement = getDisplacement((*p)->position, (*q)->position);
      switch((*p)->interaction) {
      default:
      case 0:
	hardDiskRepulsion_sym(**p, **q, displacement, n, s);
	break;
      case 1:
	LJinteraction_sym(**p, **q, displacement, n, s);
	break;
      }
    }
  }
}

void Sectorization::wallInteractions() {
  if (doWallNeighbors) {
    for (auto pr : wallNeighbors) {
      double n=0, s=0;
      Particle *p = pr.first;
      for (auto w : pr.second) {
	vect<> displacement = getDisplacement(p->position, w->left);
	switch (p->interaction) {
	default:
	case 0:
	  hardDiskRepulsion_wall(*p, *w, displacement, n, s);
	  break;
	case 1:
	  LJinteraction_wall(*p, *w, displacement, n, s);
	  break;
	}
      }
    }
  }
  else { 
    for (auto& w : walls)
      for (int i=0; i<size; ++i) {
	Particle &p = particles[i];
	double n=0, s=0;
	vect<> displacement = getDisplacement(p.position, w.left);
	switch (p.interaction) {
	default:
	case 0:
	  hardDiskRepulsion_wall(p, w, displacement, n, s);
	  break;
	case 1:
	  LJinteraction_wall(p, w, displacement, n, s);
	  break;
	}
      }
  }
}

void Sectorization::update() {
  // Reset particle's force recordings
  for (int i=0; i<size; ++i) {
    Particle &p = particles[i];
    p.force = Zero;
    p.torque = 0;
  }
  // Half-kick velocity update, update position
  double dt = 0.5 * epsilon;
  for (int i=0; i<size; ++i) {
    Particle &p = particles[i];
    double mass = 1./p.invMass;
    p.velocity += dt * p.invMass * p.force;
    p.omega    += dt * p.invII * p.torque;
    p.position += epsilon * p.velocity;    
    wrap(p.position);
    // Apply gravity (part of step 3)
    p.force += (gravity*mass - drag*p.velocity);
    // Apply temperature force (part of step 3)
    if (temperature>0 && tempDelay<time-lastTemp) {
      static double DT1 = temperature/(6*viscosity*PI);
      double DT = DT1*(1./p.sigma); // Assumes Kb = 1;
      double coeffD = sqrt(2*DT)*sqrtTempDelay;
      vect<> TForce = coeffD*(randNormal()*randV()); // Addative gaussian white noise
      p.force += TForce;
    }
  }
  // Reset last temp
  if (tempDelay<time-lastTemp) lastTemp = time;
  // Update sectorization, keep track of particles that need to migrate to other processors, send particles to the correct processor
  updateSectors();
  if (buildDelay<=itersSinceBuild) {
    itersSinceBuild = 0;
    createNeighborLists();
    if (doWallNeighbors) createWallNeighborList();
  }
  // Interaction forces (step 3)
  particleInteractions();
  wallInteractions();
  // Do behaviors (particle characteristics)
  // ---------- Behaviors ----------
  // Velocity update part two (step four)
  for (int i=0; i<size; ++i) {
    Particle &p = particles[i];
    p.velocity += dt * p.invMass * p.force;
    p.omega    += dt * p.invII * p.torque;
  }
  // Update counter
  itersSinceBuild++;
  time += epsilon;
}

void Sectorization::updateSectors() {
  if (doInteractions==false || sectors==0) return;
  // Move particles to the appropriate sectors, if they leave the domain, take note so we can migrate them. Only required to update particles in the sectors we manage (i.e. not the edge sectors)
  list<Particle*> moveUp, moveDown, moveLeft, moveRight;
  for (int y=1; y<nsy-1; ++y) 
    for (int x=1; x<nsx-1; ++x ) {
      int sec = y*nsx + x;
      vector<list<Particle*>::iterator> remove;
      for (auto p=sectors[sec].begin(); p!=sectors[sec].end(); ++p) {
	vect<> position = (*p)->position;
	if (!bounds.contains(position)) ; /*{ 
	  // The particle has migrated out of the domain. This should never happen if this is the only processor
	  remove.push_back(p);
	  if (position.x<bounds.left) moveLeft.push_back(*p);
	  else if (position.y<bounds.bottom) moveDown.push_back(*p);
	  else if (bounds.right<position.x) moveRight.push_back(*p);
	  else if (bounds.top<position.y) moveUp.push_back(*p);
	} */
	else { // The particle is still in the domain
	  int secNum = getSec(position);
	  if (secNum != sec) { // Needs to change sector
	    remove.push_back(p);
	    sectors[secNum].push_back(*p);
	  }
	  else; // Doesn't need to change sector
	}
      } // End loop over sector particles
      // Remove particles from sector as neccessary
      for (auto &p : remove) sectors[sec].erase(p);
    }
  
  // If we are the only processor, we are done
  if (numProc==1) return; 
  // Wait for everyone to finish
  // COMM_WORLD.Barrier(); //**

  // Send particles to other domains as neccessary
  /*
  if (1<numProc) {
    auto start = clock();
    passParticles( 1, 0, moveRight);
    passParticles(-1, 0, moveLeft);
    passParticles(0,  1, moveUp);
    passParticles(0, -1, moveDown);
    auto end = clock();
    transferTime += (double)(end-start)/CLOCKS_PER_SEC;
  }
  */
  // Tell other domains what is happening at boundary sectors
  //**---------------------------
  
  return;
}

void Sectorization::addParticle(Particle p) {
  plist.push_back(p);
}

void Sectorization::addWall(Wall w) {
  walls.push_back(w);
}

inline void Sectorization::wrap(vect<> &pos) {
  if (pos.x<simBounds.left)       pos.x = simBounds.right-fmod(simBounds.left-pos.x, simBounds.right-simBounds.left);
  else if (simBounds.right<pos.x) pos.x = fmod(pos.x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
  if (pos.y<simBounds.bottom)     pos.y = simBounds.top-fmod(simBounds.bottom-pos.y, simBounds.top-bounds.bottom);
  else if (simBounds.top<pos.y)   pos.y = fmod(pos.y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
}

inline int Sectorization::getSec(vect<> position) {
  double edge = cutoff + skinDepth;
  int sx = (position.x-bounds.left)/secWidth, sy = (position.y-bounds.bottom)/secHeight;

  // If inside the domain, we are done
  if (-1<sx && sx<nsx-2 && -1<sy && sy<nsy-2) return (sy+1)*nsx + (sx+1);
  
  // Handle X
  if (sx==-1) {
    // Case 1: This domain is the furthest to the right in the simulation bounds and the particle is in the right edge sector, which wraps to the first sector at the beginning of the simulation bounds
    if (bounds.right==simBounds.right) sx=nsx-1; // To the right of the
    // Case 2: Otherwise
    else sx=0;
  }
  else { // sx==nsx-2
    // Case 1: This domain is the furthest to left of the simulation bounds and the particle is in the left edge sector, which wraps to the last sector at the end of the simulation bounds
    if (bounds.left==simBounds.left) sx=0;
    // Case 2: Otherwise
    else sx = nsx-1;
  }
  
  // Handle Y
  if (sy==-1) {
    // Case 1: This domain is the furthest to the right in the simulation bounds and theparticle is in the right edge sector, which wraps to the first sector at the beginning of the simulation bounds
    if (bounds.top==simBounds.top) sy=nsy-1; // To the right of the
    // Case 2: Otherwise
    else sy=0;
  }
  else { // sy==nsy-2
    // Case 1: This domain is the furthest to left of the simulation bounds and the particle is in the left edge sector, which wraps to the last sector at the end of the simulation bounds
    if (bounds.bottom==simBounds.bottom) sy=0;
    // Case 2: Otherwise
    else sy = nsy-1;
  }
  
  return sy*nsx + sx;
}

inline void Sectorization::add(Particle *p) {
  if (sectors) {
    int sec = getSec(p->position);
    sectors[sec].push_back(p);
  }
}

inline void Sectorization::createNeighborLists() {
  if (sectors==0) return;
  if (!doInteractions) return;
  // The square of the neighbor distance
  double distance = sqr(cutoff+skinDepth);
  // Get rid of the old neighbor lists
  neighborList.clear();
  // Create new neighbor lists
  for (int y=1; y<nsy-1; ++y)
    for (int x=1; x<nsx-1; ++x) {
      for (auto p=sectors[y*nsx+x].begin(); p!=sectors[y*nsx+x].end(); ++p) {
	// Create a neighbor list for each particle, using the cells
	list<Particle*> nlist;
	nlist.push_back(*p); // You are at the head of the list
	// Create symmetric lists, so only check the required surrounding sectors ( * ) around the sector you are in ( <*> )
	// +---------+
	// | *  x  x |
	// | * <*> x |
	// | *  *  x |
	// +---------+
	
	// Check the sector you are in
	auto q = p; ++q;
	if (q!=sectors[y*nsx+x].end()) // Same sector
	  for (; q!=sectors[y*nsx+x].end(); ++q) 
	    nlist.push_back(*q);
	// Bottom left	    
	int sx = x-1, sy = y-1;
	for (auto &q : sectors[sy*nsx+sx]) {
	  vect<> r = getDisplacement((*p)->position, q->position);
	  if (sqr(r)<distance) nlist.push_back(q);
	}
	// Bottom
	sx = x;
	for (auto &q : sectors[sy*nsx+sx]) {
	  vect<> r = getDisplacement((*p)->position, q->position);
          if (sqr(r)<distance) nlist.push_back(q);
        }
	// Left
	sx = x-1; sy = y;
	for (auto &q : sectors[sy*nsx+sx]) {
	  vect<> r = getDisplacement((*p)->position, q->position);
          if (sqr(r)<distance) nlist.push_back(q);
        }
	// Top left
	sy = y+1;
	for (auto &q : sectors[sy*nsx+sx]) {
	  vect<> r = getDisplacement((*p)->position, q->position);
          if (sqr(r)<distance) nlist.push_back(q);
        }

	// Add the neighbor list to the collection if the particle has neighbors
	if (nlist.size()>1) neighborList.push_back(nlist);
      }
    }
}

inline void Sectorization::createWallNeighborList() {
  // Create Wall Neighbor list
  wallNeighbors.clear();
  for (int i=0; i<size; ++i) {
    Particle &p = particles[i];
    list<Wall*> lst;
    for (auto &w : walls) {
      
      vect<> displacement = getDisplacement(p.position, w.left);
      wallDisplacement(displacement, p, w);
      if (sqr(displacement)<sqr(1.25*p.sigma))
	lst.push_back(&w);
    }
    if (!lst.empty())
      wallNeighbors.push_back(pair<Particle*, list<Wall*> >(&p, lst));
  }
}

inline void Sectorization::migrateParticles() {
  if (numProc==1) return;
  // Look for particles that need to migrate
  list<Particle*> moveUp, moveDown, moveLeft, moveRight;
  for (int y=1; y<nsy-1; ++y)
    for (int x=1; x<nsx-1; ++x ) {
      int sec = y*nsx + x;
      vector<list<Particle*>::iterator> remove;
      for (auto p=sectors[sec].begin(); p!=sectors[sec].end(); ++p) {
        vect<> position = (*p)->position;
        if (!bounds.contains(position)) {
          // The particle has migrated out of the domain. This should never happen if this is the only processor
	  remove.push_back(p);
	  if (position.x<bounds.left) moveLeft.push_back(*p);
	  else if (position.y<bounds.bottom) moveDown.push_back(*p);
	  else if (bounds.right<position.x) moveRight.push_back(*p);
	  else if (bounds.top<position.y) moveUp.push_back(*p);
        }
      }
      // Remove particles from sector as neccessary
      for (auto &p : remove) sectors[sec].erase(p);
    }
  // Send particles to other domains as neccessary
  auto start = clock();
  passParticles( 1, 0, moveRight);
  passParticles(-1, 0, moveLeft);
  passParticles(0,  1, moveUp);
  passParticles(0, -1, moveDown);
  auto end = clock();
  transferTime += (double)(end-start)/CLOCKS_PER_SEC;
  
  // Tell other domains what is happening at boundary sectors
  //**---------------------------
  
}

inline vect<> Sectorization::getDisplacement(vect<> A, vect<> B) {
  // Get the correct (minimal) displacement vector pointing from B to A
  double X = A.x-B.x;
  double Y = A.y-B.y;
  if (wrapX) {
    double dx = (simBounds.right-simBounds.left)-fabs(X);
    if (dx<fabs(X)) X = X>0 ? -dx : dx;
  }
  if (wrapY) {
    double dy =(simBounds.top-simBounds.bottom)-fabs(Y);
    if (dy<fabs(Y)) Y = Y>0 ? -dy : dy;
  }
  return vect<>(X,Y);
}

inline void Sectorization::passParticles(const int tx, const int ty, list<Particle*> &allParticles) {
  int dx = rank % ndx, dy = rank / ndy;
  if (tx!=0) { // Transfer X
    int send = rank+sign(tx), recv = rank-sign(tx), sendX = dx+sign(tx), recvX = dx-sign(tx);
    if (dx%2==0) { // Even
      if (-1<sendX && sendX<ndx) passParticleSend(send, allParticles);
      if (-1<recvX && recvX<ndx) passParticleRecv(recv);
    }
    else { // Odd
      if (-1<recvX && recvX<ndx) passParticleRecv(recv);
      if (-1<sendX && sendX<ndx) passParticleSend(send, allParticles);
    }
  }
  
  else { // Transfer Y
    int send = rank+sign(ty)*ndx, recv = rank-sign(ty)*ndx, sendY = dy+sign(ty), recvY = dy-sign(ty);
    
    if (dy%2==0) { // Even
      if (-1<sendY && sendY<ndy) passParticleSend(send, allParticles);
      if (-1<recvY && recvY<ndy) passParticleRecv(recv);
    }
    else { // Odd
      if (-1<recvY && recvY<ndy) passParticleRecv(recv);
      if (-1<sendY && sendY<ndy) passParticleSend(send, allParticles);
    }
  }
  // COMM_WORLD.Barrier(); //** Neccessary ?
}

inline void Sectorization::passParticleSend(const int send, list<Particle*> &allParticles) {
  // Send
  int sz = allParticles.size();
  COMM_WORLD.Send(&sz, 1, MPI_INT, send, 0);
  Particle *buffer = new Particle[sz];
  int i=0;
  for (auto p : allParticles) {
    buffer[i] = *p;
    ++i;
  }
  COMM_WORLD.Send(buffer, sz, PARTICLE, send, 0);
  delete [] buffer;
}

inline void Sectorization::passParticleRecv(const int recv) {
  // Recieve
  int sz = -1;
  COMM_WORLD.Recv(&sz, 1, MPI_INT, recv, 0);
  Particle *buffer = new Particle[sz];
  COMM_WORLD.Recv(buffer, sz, PARTICLE, recv,0);
  // Add particles to sector
  for(int i=0; i<sz; ++i) addParticle(buffer[i]);
  delete [] buffer;
}

inline void Sectorization::updatePList() {
  plist.clear();
  for (int i=0; i<size; ++i) plist.push_back(particles[i]);
}
