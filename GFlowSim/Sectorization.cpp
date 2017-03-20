#include "Sectorization.h"

using MPI::COMM_WORLD;

Sectorization::Sectorization() : nsx(3), nsy(3), secWidth(0), secHeight(0), time(0), epsilon(1e-4), sqrtEpsilon(1e-4), transferTime(0), wrapX(false), wrapY(false), doInteractions(true), drag(0), gravity(Zero), temperature(0), viscosity(1.308e-3), tempDelay(5e-3), sqrtTempDelay(sqrt(5e-3)), lastTemp(0), positionTracker(0), size(0), asize(0), px(0), it(0), sectors(0), doWallNeighbors(true), remakeToUpdate(false), cutoff(0.1), skinDepth(0.05), itersSinceBuild(0), buildDelay(20), numProc(1), rank(0) {
  // Get our rank and the number of processors
  rank =    COMM_WORLD.Get_rank();
  numProc = COMM_WORLD.Get_size();  
  // Define the particle data type
  int sz = static_cast<int>(sizeof(Particle)/sizeof(MPI_DOUBLE)); // Should be 16
  MPI_Type_contiguous( sz, MPI_DOUBLE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );
  // Set all pointers to zero, create arrays
  zeroPointers();
  // Diffusion helper constant
  DT1 = temperature/(6*viscosity*PI);  
}

Sectorization::~Sectorization() {
  if (sectors)   delete [] sectors;
  if (px)        delete [] px;
  sectors = 0;
  px = 0;
}

void Sectorization::initialize() {
  // Reset data
  time = 0;
  transferTime = 0;
  lastTemp = 0;
  itersSinceBuild = 0;
  // Reset debug times
  firstHalfKick = 0;
  secondHalfKick = 0;
  updateSectorsTime = 0;
  neighborListTime = 0;
  wallNeighborListTime = 0;
  particleInteractionTime = 0;
  wallInteractionTime = 0;
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
  sectors = new list<int>[nsx*nsy+1];
  // Remake particles
  if (asize<1) asize = 4*plist.size(); //----
  if (asize<1) return; //---
  // Set particle data in arrays
  createArrays();
  setParticles();
  // Set position tracker
  for (int i=0; i<size; ++i)
    positionTracker[i] = vect<>(px[i], px[i+size]);
  // Add the particles to the proper sectors
  for (int i=0; i<size; ++i) { 
    int sec = getSec(px[i], py[i]);
    sectors[sec].push_back(i);
  }
  
  // Create the initial neighborlist
  createNeighborLists();
  if (doWallNeighbors) createWallNeighborList();
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
  for (int i=0; i<size; ++i) it[i] = inter;
}

void Sectorization::setASize(int s) {
  asize = s;
  createArrays();
}

void Sectorization::discard() {
  if (px) delete [] px; px = 0;
  if (py) delete [] py; py = 0;
  if (vx) delete [] vx; vx = 0;
  if (vy) delete [] vy; vy = 0;
  if (fx) delete [] fx; fx = 0;
  if (fy) delete [] fy; fy = 0;
  if (om) delete [] om; om = 0;
  if (tq) delete [] tq; tq = 0;
  if (sg) delete [] sg; sg = 0;
  if (im) delete [] im; im = 0;
  if (iI) delete [] iI; iI = 0;
  if (rp) delete [] rp; rp = 0;
  if (ds) delete [] ds; ds = 0;
  if (cf) delete [] cf; cf = 0;
  if (dg) delete [] dg; dg = 0;
  if (it) delete [] it; it = 0;
  if (ms) delete [] ms; ms = 0;
  size = 0;
  if (positionTracker) delete [] positionTracker;
  positionTracker = 0;
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
      vect<> displacement = getDisplacement(vect<>(px[*p],py[*p]), vect<>(px[*q],py[*q]));
      hardDiskRepulsion(pdata, *p, *q, asize, displacement, n, s); //**
      /*
      switch(static_cast<int>(it[*p])) {
      default:
      case 0:
	hardDiskRepulsion(pdata, *p, *q, asize, displacement, n, s);
	break;
      case 1:
	// LJinteraction(pdata, *p, *q, asize, displacement, n, s);
	break;
      }
      */
    }
  }
}

void Sectorization::wallInteractions() {
  if (doWallNeighbors)
    for (auto pr : wallNeighbors) {
      double n=0, s=0;
      //Particle *p = pr.first;
      int p = pr.first;
      for (auto w : pr.second) {
	vect<> displacement = getDisplacement(vect<>(px[p],py[p]), w->left);
	hardDiskRepulsion_wall(pdata, p, *w, asize, displacement, n, s); //**
	/*
	  switch (static_cast<int>(it[p])) {
	  default:
	  case 0:
	  hardDiskRepulsion_wall(pdata, p, *w, asize, displacement, n, s);
	  break;
	  case 1:
	  //LJinteraction_wall(*p, *w, displacement, n, s);
	  break;
	  }
	*/
      }
    }
  else
    for (auto& w : walls)
      for (int i=0; i<size; ++i) {
	double n=0, s=0;
	vect<> displacement = getDisplacement(px[i], py[i], w.left.x, w.left.y);
	switch (static_cast<int>(it[i])) {
	default:
	case 0:
	  hardDiskRepulsion_wall(pdata, i, w, asize, displacement, n, s);
	  break;
	case 1:
	  //LJinteraction_wall(p, w, displacement, n, s);
	  break;
	}
      }
}

void Sectorization::update() {
  // Half-kick velocity update, update position
  auto start = clock(); //--
  double dt = 0.5 * epsilon;
  bool doTemp = (temperature>0 && tempDelay<time-lastTemp);
  if (doTemp) {
#pragma vector aligned
#pragma simd
    for (int i=0; i<size; ++i) {
      // double mass = 1./im[i];
      vx[i] += dt*im[i]*fx[i];
      vy[i] += dt*im[i]*fy[i];
      om[i] += dt*iI[i]*tq[i];
      px[i] += epsilon*vx[i];
      py[i] += epsilon*vy[i];   
      wrap(px[i], py[i]);
      fx[i] = gravity.x*ms[i] - drag*vx[i];
      fy[i] = gravity.y*ms[i] - drag*vy[i];
      tq[i] = 0;
      // Temperature force
      double DT = DT1*(1./sg[i]); // Assumes Kb = 1;
      double coeffD = sqrt(2*DT)*sqrtTempDelay;
      vect<> TForce = coeffD*randNormal()*randV();
      fx[i] += TForce.x; fy[i] += TForce.y;
    }
  }
  else {
#pragma vector aligned
#pragma simd
    for (int i=0; i<size; ++i) {
      // double mass = 1./im[i];
      vx[i] += dt*im[i]*fx[i];
      vy[i] += dt*im[i]*fy[i];
      om[i] += dt*iI[i]*tq[i];
      px[i] += epsilon*vx[i];
      py[i] += epsilon*vy[i];

      wrap(px[i], py[i]);
      fx[i] = gravity.x*ms[i] - drag*vx[i];
      fy[i] = gravity.y*ms[i] - drag*vy[i];
      tq[i] = 0;
    }
  }
  auto end = clock();
  firstHalfKick += (double)(end-start)/CLOCKS_PER_SEC;
  // Reset last temp
  if (tempDelay<time-lastTemp) lastTemp = time;
  // Update sectorization, keep track of particles that need to migrate to other processors, send particles to the correct processor
  if (buildDelay<=itersSinceBuild) {
    itersSinceBuild = 0;
    createNeighborLists(); // Create the neighborhood lists
    start = clock(); //--
    if (doWallNeighbors) createWallNeighborList();
    end = clock(); //--
    wallNeighborListTime += (double)(end-start)/CLOCKS_PER_SEC; //--
  }
  
  // Interaction forces (step 3)
  start = clock(); //--
  particleInteractions();
  end = clock(); //--
  particleInteractionTime += (double)(end-start)/CLOCKS_PER_SEC; //--
  start = clock(); //--
  wallInteractions();
  end = clock(); //--
  wallInteractionTime += (double)(end-start)/CLOCKS_PER_SEC; //--
  // Velocity update part two (step four) -- second half-kick
  start = clock(); //--
#pragma vector aligned
#pragma simd
  for (int i=0; i<size; ++i) {
    vx[i] += dt*im[i]*fx[i];
    vy[i] += dt*im[i]*fy[i];
    om[i] += dt*iI[i]*tq[i];
  }
  end = clock(); //--
  secondHalfKick += (double)(end-start)/CLOCKS_PER_SEC; //--
  // Update counter
  itersSinceBuild++;
  time += epsilon;
}

void Sectorization::updateSectors() {
  if (doInteractions==false || sectors==0) return;

  if (remakeToUpdate) {
    for (int i=0; i<nsx*nsy; ++i) sectors[i].clear();
    for (int i=0; i<size; ++i) {
      if (bounds.contains(px[i], py[i])) {
	int sec_num = getSec(px[i], py[i]);
	sectors[sec_num].push_back(i);
      }
    }
  }
  else
    for (int y=1; y<nsy-1; ++y)
      for (int x=1; x<nsx-1; ++x) {
	int sec = nsx*y + x;
	for (auto P=sectors[sec].begin(); P!=sectors[sec].end(); ++P) {
	  int p = *P;
	  if (!bounds.contains(px[p], py[p])) ; // Particle not in the domain
	  else {
	    int sec_num = getSec(px[p], py[p]);
	    if (sec_num != sec) { // Changed sectors
	      P = sectors[sec].erase(P);
	      sectors[sec_num].push_back(p); // Could push back all the moved particles at the end
	    }
	    else; // Doesn't need to change sectors
	  } // End loop over sector particles
	}
      }
}

void Sectorization::addParticle(Particle p) {
  plist.push_back(p);
}

void Sectorization::addWall(Wall w) {
  walls.push_back(w);
}

inline void Sectorization::wrap(double &x, double &y) {
  if (x<simBounds.left)       x = simBounds.right-fmod(simBounds.left-x, simBounds.right-simBounds.left);
  else if (simBounds.right<x) x = fmod(x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
  if (y<simBounds.bottom)     y = simBounds.top-fmod(simBounds.bottom-y, simBounds.top-bounds.bottom);
  else if (simBounds.top<y)   y = fmod(y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
}

inline void Sectorization::wrap(vect<> &pos) {
  if (pos.x<simBounds.left)       pos.x = simBounds.right-fmod(simBounds.left-pos.x, simBounds.right-simBounds.left);
  else if (simBounds.right<pos.x) pos.x = fmod(pos.x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
  if (pos.y<simBounds.bottom)     pos.y = simBounds.top-fmod(simBounds.bottom-pos.y, simBounds.top-bounds.bottom);
  else if (simBounds.top<pos.y)   pos.y = fmod(pos.y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
}

inline int Sectorization::getSec(const vect<> &position) {
  return getSec(position.x, position.y);
}
 
inline int Sectorization::getSec(const double x, const double y) {
  double edge = cutoff + skinDepth;
  int sx = (x-bounds.left)/secWidth, sy = (y-bounds.bottom)/secHeight;
  
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

/*
inline void Sectorization::add(Particle *p) {
  if (sectors) {
    int sec = getSec(p->position);
    sectors[sec].push_back(p);
  }
}
*/

inline void Sectorization::createNeighborLists() {
  if (sectors==0) return;
  if (!doInteractions) return;
  
  // Find an upper bound on how far particles might have moved from one another
  double m0 = 0, m1 = 0;
  for (int i=0; i<size; ++i) {
    double distSqr = sqr(getDisplacement(positionTracker[i], vect<>(px[i], py[i])));
    if (distSqr>m0) m0 = distSqr;
    else if (distSqr>m1) m1 = distSqr;
  }
  double maxDiff = sqrt(m0) + sqrt(m1);
  if (maxDiff<skinDepth) return;
  
  // Update sectors
  auto start = clock(); //--
  updateSectors();
  auto end = clock(); //--
  updateSectorsTime += (double)(end-start)/CLOCKS_PER_SEC; //--
  start = clock(); //--
  // The square of the neighbor distance
  double distance = sqr(cutoff+skinDepth);
  // Get rid of the old neighbor lists
  neighborList.clear();
  // Create new neighbor lists
  for (int y=1; y<nsy-1; ++y)
    for (int x=1; x<nsx-1; ++x) {
      for (auto p=sectors[y*nsx+x].begin(); p!=sectors[y*nsx+x].end(); ++p) {
	// Create a neighbor list for each particle, using the cells
	list<int> nlist;
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
	  vect<> r = getDisplacement(px[*p], py[*p], px[q], py[q]);
	  if (sqr(r)<distance) nlist.push_back(q);
	}
	// Bottom
	sx = x;
	for (auto &q : sectors[sy*nsx+sx]) {
	  vect<> r = getDisplacement(px[*p], py[*p], px[q], py[q]);
          if (sqr(r)<distance) nlist.push_back(q);
        }
	// Left
	sx = x-1; sy = y;
	for (auto &q : sectors[sy*nsx+sx]) {
	  vect<> r = getDisplacement(px[*p], py[*p], px[q], py[q]);
          if (sqr(r)<distance) nlist.push_back(q);
        }
	// Top left
	sy = y+1;
	for (auto &q : sectors[sy*nsx+sx]) {
	  vect<> r = getDisplacement(px[*p], py[*p], px[q], py[q]);
          if (sqr(r)<distance) nlist.push_back(q);
        }

	// Add the neighbor list to the collection if the particle has neighbors
	if (nlist.size()>1) neighborList.push_back(nlist);
      }
    }
  end = clock(); //--
  neighborListTime += (double)(end-start)/CLOCKS_PER_SEC; //--
  // Update position tracker array
  for (int i=0; i<size; ++i) positionTracker[i] = vect<>(px[i], py[i]);
}

inline void Sectorization::createWallNeighborList() {
  // Create Wall Neighbor list
  wallNeighbors.clear();
  for (int i=0; i<size; ++i) {
    list<Wall*> lst;
    for (auto &w : walls) {
      vect<> displacement = getDisplacement(vect<>(px[i], py[i]), w.left);
      wallDisplacement(displacement, sg[i], w);
      if (sqr(displacement)<sqr(1.25*sg[i]))
	lst.push_back(&w);
    }
    if (!lst.empty())
      wallNeighbors.push_back(pair<int, list<Wall*> >(i, lst));
  }
}

inline void Sectorization::migrateParticles() {
  /*
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
  */
}

inline vect<> Sectorization::getDisplacement(vect<> A, vect<> B) {
  return getDisplacement(A.x, A.y, B.x, B.y);
}

inline vect<> Sectorization::getDisplacement(double ax, double ay, double bx, double by) {
  // Get the correct (minimal) displacement vector pointing from B to A
  double X = ax-bx;
  double Y = ay-by;
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

inline void Sectorization::updatePList() {
  plist.clear();
  // Create particles and push them into plist
  for (int i=0; i<size; ++i) {
    Particle p;
    p.position = vect<>(px[i], py[i]);
    p.velocity = vect<>(vx[i], vy[i]);
    p.force    = vect<>(fx[i], fy[i]);
    p.omega    = om[i];
    p.torque   = tq[i];
    p.sigma    = sg[i];
    p.invMass  = im[i];
    p.invII    = iI[i];
    p.repulsion = rp[i];
    p.dissipation = ds[i];
    p.coeff    = cf[i];
    p.drag     = dg[i];
    plist.push_back(p);
  }
}

inline void Sectorization::createArrays() {
  if (asize<1) return;
  if (px) delete [] px;
  if (py) delete [] py; 
  if (vx) delete [] vx;
  if (vy) delete [] vy;
  if (fx) delete [] fx;
  if (fy) delete [] fy;
  if (om) delete [] om;
  if (tq) delete [] tq;
  if (sg) delete [] sg;
  if (im) delete [] im;
  if (iI) delete [] iI;
  if (rp) delete [] rp;
  if (ds) delete [] ds;
  if (cf) delete [] cf;
  if (dg) delete [] dg;
  if (it) delete [] it;
  if (ms) delete [] ms;
  // Reallocate
  pdata[0]  = px = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[1]  = py = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[2]  = vx = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[3]  = vy = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[4]  = fx = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[5]  = fy = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[6]  = om = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[7]  = tq = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[8]  = sg = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[9]  = im = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[10] = iI = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[11] = rp = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[12] = ds = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[13] = cf = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[14] = dg = (double*)aligned_alloc(64, asize*sizeof(double));
  pdata[15] = it = (double*)aligned_alloc(64, asize*sizeof(double));
  ms = (double*)aligned_alloc(64, asize*sizeof(double)); // Mass array
  for (int i=0; i<asize; ++i) it[i] = -1.;
  // Set position tracker array
  if (positionTracker) delete [] positionTracker;
  positionTracker = (vect<>*)aligned_alloc(64, asize*sizeof(vect<>));
}
 
inline void Sectorization::zeroPointers() {
  positionTracker = 0;
  px = py = vx = vy = fx = fy = om = tq = sg = im = iI = rp = ds = cf = dg = it = ms = 0;
  sectors = 0;
}

inline void Sectorization::setParticles() {
  int i=0;
  for (auto p : plist) {
    px[i] = p.position.x;
    py[i] = p.position.y;
    vx[i] = p.velocity.x;
    vy[i] = p.velocity.y;
    fx[i] = p.force.x;
    fy[i] = p.force.y;
    om[i] = p.omega;
    tq[i] = p.torque;
    sg[i] = p.sigma;
    im[i] = p.invMass;
    iI[i] = p.invII;
    rp[i] = p.repulsion;
    ds[i] = p.dissipation;
    cf[i] = p.coeff;
    dg[i] = p.drag;
    it[i] = p.interaction;
    ms[i] = 1./p.invMass;  // Mass array
    ++i;
  }
  size = plist.size();
  for ( ; i<asize; ++i) it[i] = -1; // No particle stored here
 }

inline void Sectorization::atom_move() {
  // Assumes that only particles that were in boundary might need to move
  // list<int> tl, tm, tr;
  // list<int> ml,     mr;
  // list<int> bl, bm, br;
  list<int> move_lsts[9];

  for (int j=1; j<nsy-1; ++j)
    for (int i=1; i<nsx-1; ++i) {
      int n_sec = nsx*j+i;
      for (auto p : sectors[n_sec] ) {
	// Check if the particle left the domain
	if (!bounds.contains(vect<>(px[p], py[p]))) {
	  int x = 1, y = 1;
	  if (px[p]<bounds.left) x = 0;
	  else if (bounds.right<px[p]) x = 2;
	  if (py[p]<bounds.bottom) y=0;
	  else if (bounds.top<py[p]) y=2;
	  int n_lst = 3*y+x;
	  if (n_lst!=4) move_lsts[n_lst].push_back(p); // Push back the index of the particle that needs to move
	}
      }
    }
  
  // Do the actual migration
  passParticles(-1, -1, move_lsts[0]); // bl
  passParticles( 0, -1, move_lsts[1]); // bm
  passParticles(+1, -1, move_lsts[2]); // br
  passParticles(-1,  0, move_lsts[3]); // ml
  passParticles(+1,  0, move_lsts[5]); // mr
  passParticles(-1, +1, move_lsts[6]); // tl
  passParticles( 0, +1, move_lsts[7]); // tm
  passParticles(+1, +1, move_lsts[8]); // tr
}

inline void Sectorization::passParticles(int tx, int ty, const list<int> &allParticles) {
  // First, figure out the x,y coordinates of this sector
  int dx = rank % ndx, dy = rank / ndx;
  // Figure out what processor we pass to, do evens first, then odds
  int sx = dx+tx, sy = dy+ty, send = -1; // send = -1 will mean don't send
  if (-1<sx && sx<ndx && -1<sy && sy<ndy) send = ndx*sy+sx;
  int rx = dx-tx, ry = dy-ty, recv = -1; // recv = -1 will mean don't recieve
  if (-1<rx && rx<ndx && -1<ry && ry<ndy) recv = ndx*ry+rx;
  // Determine the "parity" of this processor
  bool even = ((ty==0 && dx%2==0) || (ty!=0 && dy%2==0)); // You are "even" if we are passing purely in the x-direction and you have even dx, or if we are not passing purely in the x-direction and you have even dy
  // Pass the particles. Even passes first, then odd
  if (even) { // EVEN
    if (-1<send) passParticleSend(send, allParticles);
    if (-1<recv) passParticleRecv(recv);
  }
  else {      // ODD
    if (-1<recv) passParticleRecv(recv);
    if (-1<send) passParticleSend(send, allParticles);   
  }

}

inline void Sectorization::passParticleSend(const int send, const list<int> &allParticles) {
  int sz = allParticles.size();
  COMM_WORLD.Send(&sz, 1, MPI_INT, send, 0); //** Isend
  double *buffer = new double[16*sz];
  int i=0;
  // Put particles into buffer
  for (auto p : allParticles) {
    buffer[16*i+0 ] = px[p];
    buffer[16*i+1 ] = py[p];
    buffer[16*i+2 ] = vx[p];
    buffer[16*i+3 ] = vy[p];
    buffer[16*i+4 ] = fx[p];
    buffer[16*i+5 ] = fy[p];
    buffer[16*i+6 ] = om[p];
    buffer[16*i+7 ] = tq[p];
    buffer[16*i+8 ] = sg[p];
    buffer[16*i+9 ] = im[p];
    buffer[16*i+10] = iI[p];
    buffer[16*i+11] = rp[p];
    buffer[16*i+12] = ds[p];
    buffer[16*i+13] = cf[p];
    buffer[16*i+14] = dg[p];
    buffer[16*i+15] = it[p];
    ++i;
    // Remove the particle from the particle list by setting its interaction to -1
    it[p] = -1;
  }
  COMM_WORLD.Send(buffer, sz*16*sizeof(double), MPI_DOUBLE, send, 0);
  delete [] buffer;
}

inline void Sectorization::passParticleRecv(const int recv) {
  int sz = -1;
  COMM_WORLD.Recv(&sz, 1, MPI_INT, recv, 0);
  double *buffer = new double[16*sz];
  COMM_WORLD.Recv(buffer, sz*16*sizeof(double), MPI_DOUBLE, recv, 0);
  int j=0;
  for (int i=0; i<sz; i++) {
    for ( ;it[j]!=-1; ++j); // Find the next open spot in which to put particle
    px[j] = buffer[16*i+0 ];
    py[j] = buffer[16*i+1 ];
    vx[j] = buffer[16*i+2 ];
    vy[j] = buffer[16*i+3 ];
    fx[j] = buffer[16*i+4 ];
    fy[j] = buffer[16*i+5 ];
    om[j] = buffer[16*i+6 ];
    tq[j] = buffer[16*i+7 ];
    sg[j] = buffer[16*i+8 ];
    im[j] = buffer[16*i+9 ];
    iI[j] = buffer[16*i+10];
    rp[j] = buffer[16*i+11];
    ds[j] = buffer[16*i+12];
    cf[j] = buffer[16*i+13];
    dg[j] = buffer[16*i+14];
    it[j] = buffer[16*i+15];
    ms[j] = 1./im[j]; // Mass array
    ++j;
  } 
  delete [] buffer;
}

inline void Sectorization::atom_copy() {
  
}

