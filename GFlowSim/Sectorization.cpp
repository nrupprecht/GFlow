#include "Sectorization.h"

using MPI::COMM_WORLD;

Sectorization::Sectorization() : nsx(3), nsy(3), secWidth(0), secHeight(0), epsilon(1e-4), wrapX(false), wrapY(false), sectors(0), cutoff(0.1), skinDepth(0.025), numProc(1), rank(0) {
  // Get our rank and the number of processors
  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();  
  // Define the particle data type
  MPI_Type_contiguous( 16, MPI_DOUBLE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );
}

void Sectorization::initialize(double cut) {
  if (cutoff>0) cutoff = cut;
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

void Sectorization::setSimBounds(Bounds b) {
  simBounds = b;
}

void Sectorization::interactions() {
  return; //** STUB

  for (int y=1; y<nsy-1; ++y) 
    for (int x=1 ; x<nsx-1 ; ++x) 
      for (auto &p : sectors[y*nsx+x]) {
	// Check for interactions with particles in the surrounding sectors
	
      }
}

void Sectorization::update() {
  double dt = 0.5 * epsilon;
  for (auto &p : particles) {
    double mass = 1./p.invMass;
    p.velocity += dt * p.invMass * p.force;
    p.position += epsilon * p.velocity;    
    wrap(p.position);
    // Apply gravity (part of step 3)
    p.force += gravity*mass;
  }
  // Update sectorization, keep track of particles that need to migrate to other processors, send particles to the correct processor
  updateSectors();
  // Reset particle's force recordings
  for (auto &p : particles) p.force = Zero;
  // Interaction forces
  interactions();
  // Do behaviors (particle characteristics)
  // ---------- Behaviors ----------
  // Velocity update part two (step four)
  for (auto &p :particles) p.velocity += dt * p.invMass * p.force;
}

void Sectorization::updateSectors() {
  if (sectors==0) return;
  // Move particles to the appropriate sectors, if they leave the domain, take note so we can migrate them
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
  
  if (1<numProc) {
    int data = rank, recv = -1;
    // Even Send
    if (rank%2==0) {
      if (rank+1<numProc) {
	int destination = rank+1;
	COMM_WORLD.Send( &data, 1, MPI_INT, destination, 0);
      }
    }
    else {
      int sender = rank-1;
      COMM_WORLD.Recv( &recv, 1, MPI_INT, sender, 0);
    }
    
    // Odd send
    if (rank%2) {
      int destination = (rank+1)%numProc;
      COMM_WORLD.Send( &data, 1, MPI_INT, destination, 0);
    }
    else {
      if (rank!=0 || numProc%2==0) {
	int sender = rank!=0 ? rank-1 : numProc-1;
	COMM_WORLD.Recv( &recv, 1, MPI_INT, sender, 0);
      }
    }
    // Wrap if numProc is odd
    if (numProc%2) {
      if (rank==numProc-1) {
	int destination = 0;
	COMM_WORLD.Send( &data, 1, MPI_INT, destination, 0);
      }
      else if (rank==0) {
	int sender = numProc-1;
	COMM_WORLD.Recv( &recv, 1, MPI_INT, sender, 0);
      }
    }
    // Write a little message
    cout << "Rank, recieved: " << rank << " " << recv << endl;
  }
}

void Sectorization::addParticle(Particle p) {
  int sec = getSec(p.position);
  particles.push_back(p);
  Particle *P = &(*particles.rbegin());
  if (sectors) sectors[sec].push_back(P);
}

inline void Sectorization::wrap(vect<> &pos) {
  if (pos.x<simBounds.left)       pos.x = simBounds.right-fmod(simBounds.left-pos.x, simBounds.right-simBounds.left);
  else if (simBounds.right<pos.x) pos.x = fmod(pos.x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
  if (pos.y<simBounds.bottom)     pos.y = bounds.top-fmod(simBounds.bottom-pos.y, simBounds.top-bounds.bottom);
  else if (simBounds.top<pos.y)   pos.y = fmod(pos.y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
}

inline int Sectorization::getSec(vect<> position) {
  double edge = cutoff + skinDepth;

  // THIS NEEDS TO CHANGE TO TAKE INTO ACCOUNT PARTICLES IN THE EDGE SECTORS (NOT HANDLED BY THIS DOMAIN, BUT CAN INTERACT WITH OTHER PARTICLES IN THIS DOMAIN)

  int sx = (position.x)*secWidth, sy = (position.y-bounds.bottom)*secHeight;
  return sy*nsx + sx;
}

inline void Sectorization::createNeighborLists() {
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
