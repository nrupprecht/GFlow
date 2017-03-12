#include "Sectorization.h"

using MPI::COMM_WORLD;

Sectorization::Sectorization() : nsx(3), nsy(3), secWidth(0), secHeight(0), epsilon(1e-4), wrapX(false), wrapY(false), doInteractions(true), drag(0), sectors(0), cutoff(0.1), skinDepth(0.025), itersSinceBuild(0), buildDelay(20), numProc(1), rank(0) {
  // Get our rank and the number of processors
  rank =    MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();  
  // Define the particle data type
  int size = sizeof(Particle)/sizeof(MPI_DOUBLE); // Should be 16
  MPI_Type_contiguous( size, MPI_DOUBLE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );
}

Sectorization::~Sectorization() {
  if (sectors) delete [] sectors;
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
  // Remake
  if (sectors) delete [] sectors;
  sectors = new list<Particle*>[nsx*nsy+1]; // +1 For special sector
  for (auto &p : particles) add(&p);
  
  // Create the initial neighborlist
  createNeighborLists();
  itersSinceBuild = 0;
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

void Sectorization::discard() {
  particles.clear();
  if (sectors) for (int i=0; i<nsx*nsy+1; ++i) sectors[i].clear();
  neighborList.clear();
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
  // Doing this fairly naievely
  for (auto& w : walls)
    for (auto& p : particles) {
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

void Sectorization::update() {
  // Reset particle's force recordings
  for (auto &p : particles) p.force = Zero;
  // Half-kick velocity update, update position
  double dt = 0.5 * epsilon;
  for (auto &p : particles) {
    double mass = 1./p.invMass;
    p.velocity += dt * p.invMass * p.force;
    p.position += epsilon * p.velocity;    
    wrap(p.position);
    // Apply gravity (part of step 3)
    p.force += (gravity*mass - drag*p.velocity);
  }
  // Update sectorization, keep track of particles that need to migrate to other processors, send particles to the correct processor
  updateSectors();
  if (buildDelay<=itersSinceBuild) {
    itersSinceBuild = 0;
    createNeighborLists();
  }
  // Interaction forces (step 3)
  particleInteractions();
  wallInteractions();
  // Do behaviors (particle characteristics)
  // ---------- Behaviors ----------
  // Velocity update part two (step four)
  for (auto &p :particles) p.velocity += dt * p.invMass * p.force;
  // Update counter
  itersSinceBuild++;
}

void Sectorization::updateSectors() {
  if (sectors==0) return;
  // Move particles to the appropriate sectors, if they leave the domain, take note so we can migrate them. Only required to update particles in the sectors we manage (i.e. not the edge sectors)
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
  COMM_WORLD.Barrier();

  // Send particles to other domains as neccessary
  //**--------------------------------------------
  return;
}

void Sectorization::addParticle(Particle p) {
  particles.push_back(p);
  Particle *P = &(*particles.rbegin());
  add(P);
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
