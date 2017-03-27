#include "Sectorization.h"

using MPI::COMM_WORLD;

Sectorization::Sectorization() : nsx(3), nsy(3), secWidth(0), secHeight(0), time(0), epsilon(1e-4), sqrtEpsilon(1e-4), transferTime(0), wrapX(false), wrapY(false), doInteractions(true), drag(0), gravity(0), temperature(0), viscosity(1.308e-3), tempDelay(5e-3), sqrtTempDelay(sqrt(5e-3)), lastTemp(0), size(0), array_end(0), asize(0), esize(0), earray_end(0), easize(0), redoLists(false), doWallNeighbors(true), remakeToUpdate(false), cutoff(0.1), skinDepth(0.025), itersSinceBuild(0), buildDelay(20), numProc(1), rank(0) {
  // Get our rank and the number of processors
  rank =    COMM_WORLD.Get_rank();
  numProc = COMM_WORLD.Get_size();  
  // Set all pointers to zero, create arrays
  zeroPointers();
  // Diffusion helper constant
  DT1 = temperature/(6*viscosity*PI);  
  // Set work comm as COMM_WORLD for now
  CommWork = COMM_WORLD;
}

Sectorization::~Sectorization() {
  if (sectors)   delete [] sectors;
  for (int i=0; i<15; ++i)
    if (pdata[i]) delete [] pdata[i];
  sectors = 0;
}

void Sectorization::initialize() {
  // Reset data
  time = 0;
  transferTime = 0;
  lastTemp = 0;
  itersSinceBuild = 0;
  // Find out what the cutoff should be -- The sum of the largest two interaction radii
  floatType mx=0, snd=0;
  for (auto &p : plist) {
    if (p.sigma>mx) mx = p.sigma;
    else if (p.sigma>snd) snd = p.sigma;
  }
  cutoff = mx+snd;
  // First estimate
  secWidth = secHeight = (cutoff+skinDepth);
  nsx = static_cast<int>( max(1., (bounds.right-bounds.left)/secWidth) ); 
  nsy = static_cast<int>( max(1., (bounds.top-bounds.bottom)/secHeight) );
  // Actual width and height
  secWidth = (bounds.right-bounds.left)/nsx; 
  secHeight = (bounds.top-bounds.bottom)/nsy; 
  // Add for edge sectors
  nsx += 2; nsy += 2;
  // Remake sectors
  if (sectors) delete [] sectors;
  sectors = new list<int>[nsx*nsy+1];
  // Remake particles
  if (asize<1) asize = 4*plist.size(); //----  Ad hoc
  if (asize<1) return;
  easize = asize; //---- Ad hoc
  // Set particle data in arrays
  createArrays();
  setParticles();
  // Set position tracker
  
  // Uncommenting this statement makes the program run at < half the speed it usually runs at
  for (int i=0; i<size; ++i) positionTracker[i] = vec2(px[i], py[i]);
  

  // Add the particles to the proper sectors
  for (int i=0; i<size; ++i) { 
    int sec = getSec(px[i], py[i]);
    sectors[sec].push_back(i);
  }
  
  // Create the initial neighborlist
  createNeighborLists(true);
  if (doWallNeighbors) createWallNeighborList();
}

double Sectorization::getAvePerNeighborList() {
  int count;
  for (auto &n : neighborList) count += n.size();
  return static_cast<double>(count)/neighborList.size();
}

list<Particle>& Sectorization::getParticles() {
  updatePList();
  return plist;
}

void Sectorization::setBounds(floatType l, floatType r, floatType b, floatType t) {
  bounds = Bounds(l, r, b, t);
}

void Sectorization::setBounds(Bounds b) {
  bounds = b;
}

void Sectorization::setSimBounds(floatType l, floatType r, floatType b, floatType t) {
  simBounds = Bounds(l, r, b, t);
}

void Sectorization::setInteractionType(int inter) {
  for (int i=0; i<size; ++i) it[i] = inter;
}

void Sectorization::setASize(int s) {
  asize = s;
  easize = s; //-- AD HOC
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
  if (th) delete [] th; th = 0;
  if (it) delete [] it; it = 0;
  if (ms) delete [] ms; ms = 0;
  for (int i=0; i<15; ++i) pdata[i] = 0;
  size = 0;
  array_end = 0;
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
  for (auto &nl : neighborList) {
    auto p = nl.begin();       // The particle whose list this is
    int i = *p;
    if (it[i]<0) continue;     // This particle is gone
    auto q = p; ++q;           // Start with the particle after the head particle
    int ibase = it[i]<<2;
    for (; q!=nl.end(); ++q) { // Try to interact with all other particles in the nl
      int j = *q;
      vec2 displacement = vec2(px[j]-px[i], py[j]-py[i]); //** getDisplacement(vec2(px[j],py[j]), vec2(px[i],py[i]));
      floatType Fn=0, Fs=0;
      int iType = ibase+it[j]; // Same as 4*it[i]+it[j]
      switch(iType) {
	// Both are hard disks
      case 0:
	hardDiskRepulsion(pdata, i, j, asize, displacement, Fn, Fs);
	break;
      case 1: // Sphere - LJ --> LJ
	LJinteraction(pdata, i, j, asize, displacement, Fn, Fs);
	break;
      case 2: // Sphere - Triangle
	// UNIMPLEMENTED
	break;
      case 3: // UNIMPLEMENTED
	hardDiskRepulsion(pdata, i, j, asize, displacement, Fn, Fs);
	break;
      case 4: // LJ - Sphere --> LJ
      case 5: // LJ - LJ --> LJ
      case 6: // LJ - Triangle --> LJ
      case 7: // UNIMPLEMENTED
	LJinteraction(pdata, i, j, asize, displacement, Fn, Fs);
	break;
      case 8: // Triangle - Sphere
	// UNIMPLEMENTED
	break;
      case 9: // Triangle - LJ --> LJ
	LJinteraction(pdata, i, j, asize, displacement, Fn, Fs);
	break;
      case 10: // Triangle - Triangle
	TriTriInteraction(pdata, i, j, asize, displacement, Fn, Fs);
	break;
      default:
	break;	
      } // End outer switch
    }
  }
}

void Sectorization::wallInteractions() {
  if (doWallNeighbors)
    for (auto pr : wallNeighbors) {
      int i = pr.first;
      if (it[i]<0) continue; // Particle has moved out
      for (auto w : pr.second) {
	floatType Fn=0, Fs=0;
	vec2 displacement = getDisplacement(px[i], py[i], w->left.x, w->left.y);
	switch (it[i]) {
        default:
        case 0:
          hardDiskRepulsion_wall(pdata, i, *w, asize, displacement, Fn, Fs);
          break;
        case 1:
          hardDiskRepulsion_wall(pdata, i, *w, asize, displacement, Fn, Fs);
	  //LJinteraction_wall(*p, *w, displacement, Fn, Fs);
          break;
        }
      }
    }
  else
    for (auto& w : walls)
      for (int i=0; i<array_end; ++i) {
	if (it[i]<0) continue;
	floatType Fn=0, Fs=0;
	vec2 displacement = getDisplacement(px[i], py[i], w.left.x, w.left.y);
	switch (it[i]) {
	default:
	case 0:
	  hardDiskRepulsion_wall(pdata, i, w, asize, displacement, Fn, Fs);
	  break;
	case 1:
	  // LJinteraction_wall(p, w, displacement, Fn, Fs);
	  break;
	}
      }
}

void Sectorization::update() {
  // Half-kick velocity update, update position
  firstHalfKick();
  // Reset last temp
  if (tempDelay<time-lastTemp) lastTemp = time;
  // MPI coordination
  if (doInteractions && 1<numProc) {
    atom_move();
    // atom_copy();
  }
  // Update sectorization, keep track of particles that need to migrate to other processors, send particles to the correct processor
  if (redoLists || buildDelay<=itersSinceBuild) {
    itersSinceBuild = 0;
    // Get rid of holes in the array
    compressArrays();
    // Create the neighborhood lists
    createNeighborLists(); 
    if (doWallNeighbors) createWallNeighborList();
    redoLists = false;
  }
  // Interaction forces (step 3)
  particleInteractions();
  wallInteractions();
  // Velocity update part two (step four) -- second half-kick
  secondHalfKick(); 
  // Update counter
  itersSinceBuild++;
  time += epsilon;
}

void Sectorization::updateSectors() {
  if (doInteractions==false || sectors==0) return;
  if (remakeToUpdate) { // Totally discard and remake
    for (int i=0; i<nsx*nsy; ++i) sectors[i].clear();
    for (int i=0; i<size; ++i) {
      if (bounds.contains(px[i], py[i])) {
	int sec_num = getSec(px[i], py[i]);
	sectors[sec_num].push_back(i);
      }
      else; // Particle not in the domain
    }
    for (int i=asize; i<earray_end; ++i) {
      int sec_num = getSec(px[i], py[i]);
      sectors[sec_num].push_back(i);
    }
  }
  else
    for (int y=1; y<nsy-1; ++y)
      for (int x=1; x<nsx-1; ++x) {
	int sec = nsx*y + x;
	for (auto P=sectors[sec].begin(); P!=sectors[sec].end(); ++P) {
	  int p = *P;
	  if (it[p]==-1.) P = sectors[sec].erase(P);
	  else if (bounds.contains(px[p], py[p])) {
	    int sec_num = getSec(px[p], py[p]);
	    if (sec_num != sec) { // Changed sectors
	      P = sectors[sec].erase(P);
	      sectors[sec_num].push_back(p); // Could push back all the moved particles at the end
	    }
	    else; // Doesn't need to change sectors
	  }
	  else;  // Particle not in the domain
	} // End loop over sector particles
      }
}

void Sectorization::addParticle(Particle p) {
  plist.push_back(p);
}

void Sectorization::addWall(Wall w) {
  walls.push_back(w);
}

string Sectorization::printSectors() {
  stringstream stream;
  string str;
  stream << "{";
  for (int y=0; y<nsy; ++y) {
    stream << "{";
    for (int x=0; x<nsx; ++x) {
      stream << sectors[nsx*y+x].size();
      if (x!=nsx-1) stream << ",";
    }
    stream << "}";
    if (y!=nsy-1) stream << ",";
  }
  stream << "}";
  stream >> str;
  return str;
}

inline void Sectorization::firstHalfKick() {
  double dt = 0.5 * epsilon;
  bool doTemp = (temperature>0 && tempDelay<time-lastTemp);
  if (doTemp) {
#pragma vector aligned
#pragma simd
    for (int i=0; i<array_end; ++i) {
      // floatType mass = 1./im[i];
      vx[i] += dt*im[i]*fx[i];
      vy[i] += dt*im[i]*fy[i];
      px[i] += epsilon*vx[i];
      py[i] += epsilon*vy[i];
      fx[i] = gravity.x*ms[i] - drag*vx[i];
      fy[i] = gravity.y*ms[i] - drag*vy[i];
      wrap(px[i], py[i]); // Wrap position
      om[i] += dt*iI[i]*tq[i];
      th[i] += epsilon*om[i];
      wrap(th[i]);        // Wrap theta
      tq[i] = 0;
      // Temperature force
      floatType DT = DT1*(1./sg[i]); // Assumes Kb = 1;
      floatType coeffD = sqrt(2*DT)*sqrtTempDelay;
      vec2 TForce = coeffD*randNormal()*randV();
      fx[i] += TForce.x; fy[i] += TForce.y;
    }
  }
  else {
#pragma vector aligned
#pragma simd
    for (int i=0; i<array_end; ++i) {
      // floatType mass = 1./im[i];
      vx[i] += dt*im[i]*fx[i];
      vy[i] += dt*im[i]*fy[i];
      px[i] += epsilon*vx[i];
      py[i] += epsilon*vy[i];
      fx[i] = gravity.x*ms[i] - drag*vx[i];
      fy[i] = gravity.y*ms[i] - drag*vy[i];
      wrap(px[i], py[i]);
      om[i] += dt*iI[i]*tq[i];
      th[i] += epsilon*om[i];
      tq[i] = 0;
    }
  }
}

inline void Sectorization::secondHalfKick() {
  double dt = 0.5 * epsilon;
#pragma vector aligned
#pragma simd
  for (int i=0; i<array_end; ++i) {
    vx[i] += dt*im[i]*fx[i];
    vy[i] += dt*im[i]*fy[i];
    om[i] += dt*iI[i]*tq[i];
  }
}

inline void Sectorization::wrap(floatType &x, floatType &y) {
  if (x<simBounds.left)       x = simBounds.right-fmod(simBounds.left-x, simBounds.right-simBounds.left);
  else if (simBounds.right<x) x = fmod(x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
  if (y<simBounds.bottom)     y = simBounds.top-fmod(simBounds.bottom-y, simBounds.top-bounds.bottom);
  else if (simBounds.top<y)   y = fmod(y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
}

inline void Sectorization::wrap(vec2 &pos) {
  if (pos.x<simBounds.left)       pos.x = simBounds.right-fmod(simBounds.left-pos.x, simBounds.right-simBounds.left);
  else if (simBounds.right<pos.x) pos.x = fmod(pos.x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
  if (pos.y<simBounds.bottom)     pos.y = simBounds.top-fmod(simBounds.bottom-pos.y, simBounds.top-bounds.bottom);
  else if (simBounds.top<pos.y)   pos.y = fmod(pos.y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
}

inline void Sectorization::wrap(floatType &theta) {
  if (theta<0) theta = 2*PI-fmod(-theta, 2*PI);
  else if (2*PI<theta) theta = fmod(theta, 2*PI);
}

inline int Sectorization::getSec(const vec2 &position) {
  return getSec(position.x, position.y);
}
 
inline int Sectorization::getSec(const floatType x, const floatType y) {
  floatType edge = cutoff + skinDepth;
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

inline void Sectorization::createNeighborLists(bool force) {
  if (sectors==0) return;
  if (!doInteractions) return;
  
  // Find an upper bound on how far particles might have moved from one another, unless we are forcing a redo
  if (!force) {
    floatType m0 = 0, m1 = 0;
    for (int i=0; i<array_end; ++i) {
      if (0<=it[i]) {
	floatType distSqr = sqr(getDisplacement(positionTracker[i], vec2(px[i], py[i])));
	if (distSqr>m0) m0 = distSqr;
	else if (distSqr>m1) m1 = distSqr;
      }
    }
    maxNLDiff = sqrt(m0) + sqrt(m1);
    if (maxNLDiff<skinDepth) return;
  }
  // Update sectors
  updateSectors();
  // The square of the neighbor distance
  floatType distance = sqr(cutoff+skinDepth);
  // Get rid of the old neighbor lists
  neighborList.clear();
  // Create new neighbor lists
  for (int y=1; y<nsy-1; ++y)
    for (int x=1; x<nsx-1; ++x) {
      for (auto p=sectors[y*nsx+x].begin(); p!=sectors[y*nsx+x].end(); ++p) {
	// Only treat valid particles
	int i = *p;
	if (it[i]<0) continue;
	// Create a neighbor list for each particle, using the cells
	list<int> nlist;
	nlist.push_back(i); // You are at the head of the list
	floatType sigma = sg[i], velocity = sqrt(sqr(vx[i])+sqr(vy[i]));
	floatType range = sigma * max(velocity, 1.);
	// Create symmetric lists, so only check the required surrounding sectors ( * ) around the sector you are in ( <*> )
	// +---------+
	// | *  x  x |
	// | * <*> x |
	// | *  *  x |
	// +---------+
	
	// Check the sector you are in
	auto q = p; ++q;
        if (q!=sectors[y*nsx+x].end()) // Same sector
          for (; q!=sectors[y*nsx+x].end(); ++q) {
            vec2 r = getDisplacement(px[i], py[i], px[*q], py[*q]);
            if (sqr(r)<sqr(range+sg[*q]+skinDepth)) nlist.push_back(*q);
          }
        // Bottom left
        int sx = x-1, sy = y-1;
	if (0<sy) {
	  if (0<sx) // Don't hit edge sectors
	    for (auto &q : sectors[sy*nsx+sx]) {
	      vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	      if (sqr(r)<sqr(range+sg[q]+skinDepth)) nlist.push_back(q);
	    }
	  // Bottom
	  sx = x;
	  for (auto &q : sectors[sy*nsx+sx]) {
	    vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	    if (sqr(r)<sqr(range+sg[q]+skinDepth)) nlist.push_back(q);
	  }
	}
        // Left
        sx = x-1; sy = y;
	if (0<sx) {
	  for (auto &q : sectors[sy*nsx+sx]) {
	    vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	    if (sqr(r)<sqr(range+sg[q]+skinDepth)) nlist.push_back(q);
	  }
	  // Top left
	  sy = y+1;
	  if (sy<nsy-1)
	    for (auto &q : sectors[sy*nsx+sx]) {
	      vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	      if (sqr(r)<sqr(range+sg[q]+skinDepth)) nlist.push_back(q);
	    }
	}
	// Add the neighbor list to the collection if the particle has neighbors
	if (nlist.size()>1) neighborList.push_back(nlist);
      }
    }
  // Update position tracker array
  for (int i=0; i<array_end; ++i) positionTracker[i] = vec2(px[i], py[i]);
}

inline void Sectorization::createWallNeighborList() {
  // Create Wall Neighbor list
  wallNeighbors.clear();
  for (int i=0; i<array_end; ++i) {
    if (it[i]<0) continue;
    list<Wall*> lst;
    for (auto &w : walls) {
      vec2 displacement = getDisplacement(vec2(px[i], py[i]), w.left);
      wallDisplacement(displacement, sg[i], w);
      if (sqr(displacement)<sqr(1.25*sg[i]))
	lst.push_back(&w);
    }
    if (!lst.empty())
      wallNeighbors.push_back(pair<int, list<Wall*> >(i, lst));
  }
}

inline vec2 Sectorization::getDisplacement(vec2 A, vec2 B) {
  return getDisplacement(A.x, A.y, B.x, B.y);
}

inline vec2 Sectorization::getDisplacement(floatType ax, floatType ay, floatType bx, floatType by) {
  // Get the correct (minimal) displacement vector pointing from B to A
  floatType X = ax-bx;
  floatType Y = ay-by;
  if (wrapX) {
    floatType dx = (simBounds.right-simBounds.left)-fabs(X);
    if (dx<fabs(X)) X = X>0 ? -dx : dx;
  }
  if (wrapY) {
    floatType dy =(simBounds.top-simBounds.bottom)-fabs(Y);
    if (dy<fabs(Y)) Y = Y>0 ? -dy : dy;
  }
  return vec2(X,Y);
}

inline void Sectorization::updatePList() {
  plist.clear();
  // Create particles and push them into plist
  for (int i=0; i<array_end; ++i) {
    if (it[i]<0) continue;
    Particle p;
    // Set particle data
    p.position = vec2(px[i], py[i]);
    p.velocity = vec2(vx[i], vy[i]);
    p.force    = vec2(fx[i], fy[i]);
    p.theta    = th[i];
    p.omega    = om[i];
    p.torque   = tq[i];
    p.sigma    = sg[i];
    p.invMass  = im[i];
    p.invII    = iI[i];
    p.repulsion = rp[i];
    p.dissipation = ds[i];
    p.coeff    = cf[i];
    p.interaction = it[i];
    // Add particle to list
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
  if (th) delete [] th;
  if (om) delete [] om;
  if (tq) delete [] tq;
  if (sg) delete [] sg;
  if (im) delete [] im;
  if (iI) delete [] iI;
  if (rp) delete [] rp;
  if (ds) delete [] ds;
  if (cf) delete [] cf;
  if (it) delete [] it;
  if (ms) delete [] ms;
  // Reallocate
  int tsize = asize + easize;
  pdata[0]  = px = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[1]  = py = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[2]  = vx = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[3]  = vy = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[4]  = fx = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[5]  = fy = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[6]  = th = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[7]  = om = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[8]  = tq = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[9]  = sg = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[10] = im = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[11] = iI = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[12] = rp = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[13] = ds = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  pdata[14] = cf = (floatType*)aligned_alloc(64, tsize*sizeof(floatType));
  it             =    (int*)   aligned_alloc(64, tsize*sizeof(int));
  memset(it, -1, tsize*sizeof(int)); // Each int is 4 bytes
  ms = (floatType*)aligned_alloc(64, asize*sizeof(floatType));
  // Set position tracker array
  if (positionTracker) delete [] positionTracker;
  positionTracker = (vec2*)aligned_alloc(64, asize*sizeof(vec2));
  size = 0; esize = 0;
  array_end = 0; earray_end = asize;
}
 
inline void Sectorization::zeroPointers() {
  positionTracker = 0;
  px = py = vx = vy = fx = fy = th = om = tq = sg = im = iI = rp = ds = cf = ms = 0;
  it = 0;
  for (int i=0; i<15; ++i) pdata[i] = 0;
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
    th[i] = p.theta;
    om[i] = p.omega;
    tq[i] = p.torque;
    sg[i] = p.sigma;
    im[i] = p.invMass;
    iI[i] = p.invII;
    rp[i] = p.repulsion;
    ds[i] = p.dissipation;
    cf[i] = p.coeff;
    it[i] = p.interaction;
    ms[i] = 1./p.invMass;  // Mass array
    ++i;
  }
  size = i;
  array_end = i;
  for ( ; i<asize; ++i) it[i] = -1; // No particle stored here
 }

inline void Sectorization::atom_move() {
  // Assumes that only particles that were in boundary might need to move
  //  tl, tm, tr;
  //  ml, **  mr;
  //  bl, bm, br;
  list<int> move_lsts[9]; // Entry 4 will always be unused
  // vector<int> holes;      // Holes that open because of particles leaving. We will fill these holes

  for (int i=0; i<array_end; ++i) {
    if (it[i]<0) continue;
    // Check if the particle left the domain
    if (!bounds.contains(vec2(px[i], py[i]))) {
      holes.push_back(i);
      int x = 1, y = 1;
      if (px[i]<bounds.left) x = 0;
      else if (bounds.right<px[i]) x = 2;
      if (py[i]<bounds.bottom) y=0;
      else if (bounds.top<py[i]) y=2;
      int n_lst = 3*y+x;
      if (n_lst!=4) move_lsts[n_lst].push_back(i); // Push back the index of the particle that needs to move
    }
  }

  /*
  for (int j=1; j<nsy-1; ++j)
    for (int i=1; i<nsx-1; ++i) {
      int n_sec = nsx*j+i;
      for (auto P=sectors[n_sec].begin(); P!=sectors[n_sec].end(); ++P) {
	int p = *P;
	// Check if the particle left the domain
	if (!bounds.contains(vec2(px[p], py[p]))) {
	  P = sectors[n_sec].erase(P);
	  holes.push_back(p);
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
  */
  
  // Do the actual migration
  auto start = clock();
  passParticles(-1, -1, move_lsts[0]); // bl
  passParticles( 0, -1, move_lsts[1]); // bm
  passParticles(+1, -1, move_lsts[2]); // br
  passParticles(-1,  0, move_lsts[3]); // ml
  passParticles(+1,  0, move_lsts[5]); // mr
  passParticles(-1, +1, move_lsts[6]); // tl
  passParticles( 0, +1, move_lsts[7]); // tm
  passParticles(+1, +1, move_lsts[8]); // tr
  auto end = clock();
  transferTime += (double)(end-start)/CLOCKS_PER_SEC;

}

inline void Sectorization::passParticles(int tx, int ty, const list<int> &allParticles, bool edgeParticles) {
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
    if (-1<send) passParticleSend(send, allParticles, edgeParticles);
    if (-1<recv) passParticleRecv(recv, edgeParticles);
  }
  else {      // ODD
    if (-1<recv) passParticleRecv(recv, edgeParticles);
    if (-1<send) passParticleSend(send, allParticles, edgeParticles);   
  }

}

inline void Sectorization::passParticleSend(const int send, const list<int> &allParticles, bool noErase) {
  // Send expected size
  int sz = allParticles.size();
  CommWork.Send(&sz, 1, MPI_INT, send, 0); //** Isend
  // If there are particles to send
  if (0<sz) {
    double *buffer = new double[16*sz];
    int i=0;
    // Put particles into buffer
    for (auto j : allParticles) {
      buffer[i+0 ] = px[j];
      buffer[i+1 ] = py[j];
      buffer[i+2 ] = vx[j];
      buffer[i+3 ] = vy[j];
      buffer[i+4 ] = fx[j];
      buffer[i+5 ] = fy[j];
      buffer[i+6 ] = th[j];
      buffer[i+7 ] = om[j];
      buffer[i+8 ] = tq[j];
      buffer[i+9 ] = sg[j];
      buffer[i+10] = im[j];
      buffer[i+11] = iI[j];
      buffer[i+12] = rp[j];
      buffer[i+13] = ds[j];
      buffer[i+14] = cf[j];
      buffer[i+15] = static_cast<floatType>(it[j]);
      i+=16;
      // Remove the particle from the particle list by setting its interaction to -1
      if (!noErase) it[j] = -1;
    }
    size -= sz; // Adjust our size
    // Send our data
    int msz = sz*16;
    CommWork.Send(buffer, msz, MPI_DOUBLE, send, 0); // ISend
    delete [] buffer;
  }
}

inline void Sectorization::passParticleRecv(const int recv, bool edgeParticles) {
  // Recieve expected size
  int sz = -1;
  CommWork.Recv(&sz, 1, MPI_INT, recv, 0);
  // If there are particles to recieve
  if (0<sz) {
    // Get our data
    double *buffer = new double[16*sz]; // (double*)aligned_alloc(64, 16*sz*sizeof(double));
    int msz = sz*16;
    MPI::Status status;
    CommWork.Recv(buffer, msz, MPI_DOUBLE, recv, 0, status);
    // Add particles
    int &end = edgeParticles ? earray_end : array_end;
    for (int i=0; i<sz; i++) {
      px[end] = buffer[16*i+0 ];
      py[end] = buffer[16*i+1 ];
      vx[end] = buffer[16*i+2 ];
      vy[end] = buffer[16*i+3 ];
      fx[end] = buffer[16*i+4 ];
      fy[end] = buffer[16*i+5 ];
      th[end] = buffer[16*i+6 ];
      om[end] = buffer[16*i+7 ];
      tq[end] = buffer[16*i+8 ];
      sg[end] = buffer[16*i+9 ];
      im[end] = buffer[16*i+10];
      iI[end] = buffer[16*i+11];
      rp[end] = buffer[16*i+12];
      ds[end] = buffer[16*i+13];
      cf[end] = buffer[16*i+14];
      it[end] = static_cast<int>(buffer[16*i+15]);
      ms[end] = 1./im[end]; // Mass array
      // Add to sectors -- For some reason this creates errors when using  multiple processors
      int num_sec = getSec(px[end], py[end]);
      sectors[num_sec].push_back(end);
      // Increment array counter
      ++end;
    } 
    // Adjust our size
    if (edgeParticles) esize += sz;
    else size += sz;
    // Free memory
    delete [] buffer;
    // We should redo lists since we have recieved particles
    redoLists = true;
  }
}

inline void Sectorization::compressArrays() {
  // Holes will not in general be in order
  if (holes.empty()) return;
  for (auto J=holes.begin(); J!=holes.end(); ++J) {
    int j = *J;
    while (it[array_end-1]<0) array_end--;
    if (array_end-1<=j) continue; 
    px[j] = px[array_end-1];
    py[j] = py[array_end-1];
    vx[j] = vx[array_end-1];
    vy[j] = vy[array_end-1];
    fx[j] = fx[array_end-1];
    fy[j] = fy[array_end-1];
    th[j] = th[array_end-1];
    om[j] = om[array_end-1];
    tq[j] = tq[array_end-1];
    sg[j] = sg[array_end-1];
    im[j] = im[array_end-1];
    iI[j] = iI[array_end-1];
    rp[j] = rp[array_end-1];
    ds[j] = ds[array_end-1];
    cf[j] = cf[array_end-1];
    it[j] = it[array_end-1];
    ms[j] = ms[array_end-1]; // Mass array
    it[array_end-1] = -1;    // This entry is now empty
    // Increment array end
    --array_end;
  }
  holes.clear();
}

inline void Sectorization::atom_copy() {
  // Clear edge sectors
  for (int x=0; x<nsx; ++x)   sectors[x].clear(); // Bottom edge
  for (int x=0; x<nsx; ++x)   sectors[nsx*nsy-x-1].clear(); // Top edge
  for (int y=1; y<nsy-1; ++y) sectors[nsx*y].clear(); // Left edge
  for (int y=1; y<nsy-1; ++y) sectors[nsx*y+nsx-1].clear(); // Right edge
  earray_end = asize; // Reset edge particles
  esize = 0;
  // Domain coordinates
  int dx = rank % ndx, dy = rank / ndx;
  // Arrays
  list<int> leftEdge, rightEdge, topEdge, bottomEdge;
  // Start timing
  auto start = clock();
  // Pass left edge to the left
  if (0<dx) // Or wrap y
    for (int y=1; y<nsy-1; ++y)
      for (auto i : sectors[nsx*y+1]) 
	if (0<=it[i]) leftEdge.push_back(i);
  // passParticles(-1, 0, leftEdge, true);
  passParticles(-1, 0, leftEdge);

  // Pass right edge to the right
  if (dx<ndx-1) // Or wrap x
    for (int y=1; y<nsy-1; ++y)
      for (auto i : sectors[nsx*y+nsx-1]) 
	if (0<=it[i]) leftEdge.push_back(i);
  // passParticles(+1, 0, rightEdge, true);
  passParticles(+1, 0, rightEdge);

  // Pass top edge upwards
  if (dy<ndy-1)
    for (int x=0; x<nsx; ++x)
      for (auto i : sectors[nsx*(nsy-1)+x]) 
	if (0<=it[i]) topEdge.push_back(i);
  // passParticles(0, +1, topEdge, true);
  cout << "Rank: " << rank << ", " << topEdge.size() << endl; //**
  passParticles(0, +1, topEdge);

  // Pass bottom edge downwards
  if (0<dy)
    for (int x=0; x<nsx; ++x)
      for(auto i: sectors[nsx+x]) 
	if (0<=it[i]) bottomEdge.push_back(i);
  // passParticles(0, -1, bottomEdge, true);
  // cout << bottomEdge.size() << endl; //**
  cout << "Rank: " << rank << ", " << bottomEdge.size() << endl; //**
  passParticles(0, -1, bottomEdge);
  // Timing
  auto end = clock();
  transferTime += (double)(end-start)/CLOCKS_PER_SEC;
}

