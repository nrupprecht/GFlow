#include "Sectorization.h"

using MPI::COMM_WORLD;

Sectorization::Sectorization() : nsx(3), nsy(3), secWidth(0), secHeight(0), time(0), epsilon(1e-4), sqrtEpsilon(1e-4), transferTime(0), wrapX(false), wrapY(false), doInteractions(true), drag(false), gravity(0), temperature(0), viscosity(1.308e-3), tempDelay(5e-3), sqrtTempDelay(sqrt(5e-3)), lastTemp(0), size(0), array_end(0), asize(0), esize(0), earray_end(0), easize(0), useCharacteristics(false), redoLists(false), doWallNeighbors(true), remakeToUpdate(false), cutoff(0.1), skinDepth(0.025), itersSinceBuild(0), buildDelay(20), numProc(1), rank(0), maxCutR(0), secCutR(0), maxNLDiff(0) {
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
  if (it) delete [] it;
  if (ch) {
    for (int i=0; i<asize; ++i)
      if (ch[i]) delete ch[i];
    delete [] ch;
  }
  sectors = 0;
}

void Sectorization::initialize() {
  // Reset data
  time = 0;
  transferTime = 0;
  lastTemp = 0;
  itersSinceBuild = 0;
  // Find out what the cutoff should be -- The sum of the largest two interaction radii
  for (auto &p : plist) {
    double r = particle_cutoff(p.sigma, p.interaction);
    if (r>maxCutR) maxCutR = r;
    else if (r>secCutR) secCutR = r;
  }
  makeSectors();
  // Remake particle arrays
  if (asize<1) asize = 2*plist.size(); //----  Ad hoc
  if (asize<1) return;
  easize = asize; //---- Ad hoc
  // Set particle data in arrays
  createArrays();
  setParticles();
  // Set position tracker
  for (int i=0; i<size; ++i) positionTracker[i] = vec2(px[i], py[i]);  
  // Add the particles to the proper sectors
  for (int i=0; i<size; ++i) { 
    int sec = getSec(px[i], py[i]);
    sectors[sec].push_back(i);
  }
  // Reallocate in a memory friendly way (using the sectorized data)
  memory_rearrange();
  
  // Create the initial neighborlist
  createNeighborLists(true);
  if (doWallNeighbors) createWallNeighborList();
}

double Sectorization::getAvePerNeighborList() {
  int count = 0, size = neighborList.size();
  for (auto &n : neighborList) count += n.size();
  return size>0 ? static_cast<double>(count)/size : 0;
}

vector<Particle>& Sectorization::getParticles() {
  updatePList();
  return plist;
}

void Sectorization::setCoeff(double c) {
  if (cf!=0)
    for (int i=0; i<array_end; ++i) cf[i] = c;
  for (auto &p : plist) p.coeff = c;
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

void Sectorization::setUseCharacteristics(bool use) {
  useCharacteristics = use;
  if (use) {
    // Delete old characteristics
    if (ch) {
      for (int i=0; i<asize; ++i)
	if (ch[i]) {
          delete ch[i];
          ch[i] = 0;
        }
      free(ch);
    }
    ch = (Characteristic**)aligned_alloc(64, asize*sizeof(Characteristic*));
    memset(ch, 0, asize*sizeof(Characteristic*));
  }
  else {
    // Clean up, delete old characteristics
    if (ch)
      for (int i=0; i<asize; ++i)
        if (ch[i]) {
          delete ch[i];
          ch[i] =0;
        }
    free(ch);
    ch = 0;
  }
}

void Sectorization::setInteractionType(int inter) {
  for (int i=0; i<size; ++i) it[i] = inter;
}

void Sectorization::setASize(int s) {
  updatePList();
  discard();
  asize = s;
  easize = s; //-- AD HOC
  createArrays();
  setParticles();
  //** Need to redo all lists
}

void Sectorization::stopParticles() {
  for (int i=0; i<array_end; ++i) vx[i] = vy[i] = 0;
}

void Sectorization::discard() {
  // Discard arrays only. Not plist or clist. Also, do not delete individual characteristics ( ch[i] )
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
  if (ch) delete [] ch; ch = 0;
  for (int i=0; i<15; ++i) pdata[i] = 0;
  size = 0;
  array_end = 0;
  if (positionTracker) delete [] positionTracker;
  positionTracker = 0;
  if (sectors) for (int i=0; i<nsx*nsy+1; ++i) sectors[i].clear();
  neighborList.clear();
  wallNeighbors.clear();
}

void Sectorization::discardAll() {
  // Discard the characteristics
  if (ch) {
    for (int i=0; i<array_end; ++i) 
      if (ch[i]) {
	delete ch[i];
	ch[i] = 0;
      }
    delete [] ch;
    ch = 0;
  }
  walls.clear();
  discard();
}

void Sectorization::setSimBounds(Bounds b) {
  simBounds = b;
}

inline void Sectorization::particleInteractions() {
  if (!doInteractions) return;
  // Neighbor list
  for (auto &nl : neighborList) {
    auto p = nl.begin();       // The particle whose list this is
    int i = *p;
    if (it[i]<0) continue;     // This particle is gone
    auto q = p; ++q;           // Start with the particle after the head particle
    int ibase = it[i]<<2;      // Represents the interaction type of particle p
    for (; q!=nl.end(); ++q) { // Try to interact with all other particles in the nl
      int j = *q;
      if (it[j]<0) continue;   // This particle is gone
      double Fn=0, Fs=0;
      interactionHelper(i, j, ibase, Fn, Fs);
    }
  }
  // Edge neighbor list
  for (auto &nl : edgeNeighborList) {
    auto p = nl.begin();       // The particle whose list this is
    int i = *p;
    if (it[i]<0) continue;     // This particle is gone
    auto q = p; ++q;           // Start with the particle after the head particle
    int ibase = it[i]<<2;      // Represents the interaction type of particle p
    for (; q!=nl.end(); ++q) { // Try to interact with all other particles in the nl
      int j = *q;
      if (it[j]<0) continue;   // This particle is gone
      double Fn=0, Fs=0;
      interactionHelper(i, j, ibase, Fn, Fs);
    }
  }
}

inline void Sectorization::interactionHelper(int i, int j, int ibase, double &Fn, double &Fs, bool update) {
  vec2 displacement = getDisplacement(vec2(px[j],py[j]), vec2(px[i],py[i]));
  int iType = ibase+it[j]; // Same as 4*it[i]+it[j]
  // Switch on interaction type
  switch(iType) {
  case 0: // Sphere - Sphere
    hardDiskRepulsion(pdata, i, j, asize, displacement, Fn, Fs, update);
    break;
  case 1: // Sphere - LJ --> LJ
    LJinteraction(pdata, i, j, asize, displacement, Fn, Fs, update);
    break;
  case 2: // Sphere - Triangle
    // UNIMPLEMENTED
    break;
  case 3: // Sphere - Inv R --> Inv R
    invRRepulsion(pdata, i, j, asize, displacement, Fn, Fs, update);
    break;
  case 4: // LJ - Sphere --> LJ
  case 5: // LJ - LJ --> LJ
  case 6: // LJ - Triangle --> LJ
  case 7: // UNIMPLEMENTED
    LJinteraction(pdata, i, j, asize, displacement, Fn, Fs, update);
    break;
  case 8: // Triangle - Sphere
    // UNIMPLEMENTED
    break;
  case 9: // Triangle - LJ --> LJ
    LJinteraction(pdata, i, j, asize, displacement, Fn, Fs, update);
    break;
  case 10: // Triangle - Triangle
    TriTriInteraction(pdata, i, j, asize, displacement, Fn, Fs, update);
    break;
  case 15: // Inv R - Inv R
    invRRepulsion(pdata, i, j, asize, displacement, Fn, Fs, update);
    break;
  default:
    break;
  } // End outer switch
}

inline void Sectorization::forceChoice(int choice, const double Fn, const double Fs, double& F) {
  // 0 -- Total force/pressure
  // 1 -- Normal force/pressure
  // 2 -- Shear force/pressure
  // 3 -- Inward force/pressure
  // 4 -- Outward force/pressure
  switch (choice) {
  default:
  case 0:
    F = Fn+Fs;
    break;
  case 1:
    F = Fn;
    break;
  case 2:
    F = Fs;
    break;
  case 3:
    F = clamp(Fn+Fs);
    break;
  case 4:
    F = clamp(-Fn-Fs);
    break;
  }
}

inline void Sectorization::wallInteractions() {
  if (doWallNeighbors) {
    for (auto pr : wallNeighbors) {
      int i = pr.first;
      if (it[i]<0) continue; // Particle has moved out
      for (auto w : pr.second) {
	double Fn=0, Fs=0;
	vec2 displacement = getDisplacement(px[i], py[i], w->left.x, w->left.y);
	switch (it[i]) {
        default:
        case 0:
          hardDiskRepulsion_wall(pdata, i, *w, asize, displacement, Fn, Fs);
          break;
        case 1:
          hardDiskRepulsion_wall(pdata, i, *w, asize, displacement, Fn, Fs); //**
	  //LJinteraction_wall(*p, *w, displacement, Fn, Fs);
          break;
        }
      }
    }

  }
  else
    for (auto& w : walls)
      for (int i=0; i<array_end; ++i) {
	if (it[i]<0) continue;
	double Fn=0, Fs=0;
	vec2 displacement = getDisplacement(px[i], py[i], w.left.x, w.left.y);
	switch (it[i]) {
	default:
	case 0:
	  hardDiskRepulsion_wall(pdata, i, w, asize, displacement, Fn, Fs);
	  break;
	case 1:
	  hardDiskRepulsion_wall(pdata, i, w, asize, displacement, Fn, Fs); //**
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
  // Characteristics
  if (ch)
    for (int i=0; i<array_end; ++i) 
      if (ch[i]) ch[i]->modify(pdata, this, i);
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
  // Compress array
  else compressArrays();
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

void Sectorization::insertParticle(Particle p) {
  // Check if we need to change our sector size
  double r = particle_cutoff(p.sigma, p.interaction);
  bool redo = false;
  if (r>maxCutR) {
    secCutR = maxCutR;
    maxCutR = r;
    redo = true;
  }
  else if (r>secCutR) {
    secCutR = r;
    redo = true;
  }
  // If we need to redo the sector sizing
  if (redo) {
    updatePList();
    plist.push_back(p);
    makeSectors();
    // Add the particles (except the new one - we add it below) to the proper sectors
    for (int i=0; i<size; ++i) {
      int sec = getSec(px[i], py[i]);
      sectors[sec].push_back(i);
    }
  }
  // Add particle to the arrays etc.
  if (asize<=array_end) throw OutOfBounds(); // To many particles - we should (well, could) actually resize everything in this case
  px[array_end] = p.position.x;
  py[array_end] = p.position.y;
  vx[array_end] = p.velocity.x;
  vy[array_end] = p.velocity.y;
  fx[array_end] = p.force.x;
  fy[array_end] = p.force.y;
  th[array_end] = p.theta;
  om[array_end] = p.omega;
  tq[array_end] = p.torque;
  sg[array_end] = p.sigma;
  im[array_end] = p.invMass;
  iI[array_end] = p.invII;
  rp[array_end] = p.repulsion;
  ds[array_end] = p.dissipation;
  cf[array_end] = p.coeff;
  it[array_end] = p.interaction;
  // Mark in position tracker
  positionTracker[array_end] = p.position;
  // Add the new particle to the sectorization
  int num_sec = getSec(p.position);
  sectors[num_sec].push_back(array_end);
  // Create (potentially) wall neighbor lists
  list<Wall*> lst;
  for (auto &w : walls) {
    vec2 displacement = getDisplacement(vec2(px[array_end], py[array_end]), w.left);
    wallDisplacement(displacement, sg[array_end], w);
    if (sqr(displacement)<sqr(1.25*particle_cutoff(sg[array_end], it[array_end])))
      lst.push_back(&w);
  }
  if (!lst.empty()) wallNeighbors.push_back(pair<int, list<Wall*> >(array_end, lst));
  //** Should redo neighbor lists now -- can just create a NL for the inserted particle
  
  // Increment array end
  ++array_end;
  ++size;
}

void Sectorization::insertParticle(Particle p, Characteristic *c) {
  if (array_end>=asize) throw OutOfBounds();
  if (ch) ch[array_end] = c;
  insertParticle(p);
}

void Sectorization::removeAt(int i) {
  if (i<0 || asize-1<i) throw OutOfBounds(); 
  it[i] = -1;
  // If the particle had a characteristic, destroy it
  if (ch && ch[i]) {
    delete ch[i];
    ch[i] = 0;
  }
  // There is now a hole to be filled
  holes.push_back(i);
}

void Sectorization::addWall(Wall w) {
  walls.push_back(w);
}

void Sectorization::setCharacteristic(Characteristic *C) {
  setUseCharacteristics(true);
  for (int i=0; i<array_end; ++i)
    if (it[i]!=-1) ch[i] = C->create();
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

pair<double, int> Sectorization::doStatFunction(StatFunc func) {
  // You are responsible for updating plist beforehand
  int count;
  double data = func(plist, count);
  return pair<double, int>(data, count);
}

// type : 0 - particles and walls, 1 - particles, 2 - walls
vector<PData> Sectorization::forceData(int choice, int type, bool pressure) {
  if (!doInteractions) return vector<PData>();
  // Get rid of holes in the array
  compressArrays();
  // Create the neighborhood lists
  createNeighborLists();
  // Data array with an entry for every particle
  vector<PData> data(size); // There should be no holes in the array
  // Set positions and radii
  for (int i=0; i<size; ++i) {
    std::get<0>(data.at(i)) = vec2(px[i], py[i]);
    std::get<1>(data.at(i)) = sg[i];
    std::get<2>(data.at(i)) = th[i];
    std::get<3>(data.at(i)) = it[i];
    std::get<4>(data.at(i)) = 0;
  }
  // Particle forces
  if (type==0 || type==1)
    for (auto &nl : neighborList) {
      auto p = nl.begin();       // The particle whose list this is
      int i = *p;
      if (it[i]<0) continue;     // This particle is gone
      auto q = p; ++q;           // Start with the particle after the head particle
      int ibase = it[i]<<2;      // Represents the interaction type of particle p
      for (; q!=nl.end(); ++q) { // Try to interact with all other particles in the nl
	int j = *q;
	double Fn=0, Fs=0, F=0;
	interactionHelper(i, j, ibase, Fn, Fs, false);
	forceChoice(choice, Fn, Fs, F);
	// Store
	std::get<4>(data.at(i)) += fabs(F);
	std::get<4>(data.at(j)) += fabs(F);
      }
    }
  // Wall interactions
  if (type==0 || type==2)
    for (auto& w : walls)
      for (int i=0; i<array_end; ++i) {
	if (it[i]<0) continue;
	double Fn=0, Fs=0, F=0;
	vec2 displacement = getDisplacement(px[i], py[i], w.left.x, w.left.y);
	switch (it[i]) {
	default:
	case 0:
	  hardDiskRepulsion_wall(pdata, i, w, asize, displacement, Fn, Fs);
	  break;
	case 1:
	  hardDiskRepulsion_wall(pdata, i, w, asize, displacement, Fn, Fs); //**
	  // LJinteraction_wall(p, w, displacement, Fn, Fs);
	  break;
	}
	// Add to data
	forceChoice(choice, Fn, Fs, F);    	
	std::get<4>(data.at(i)) += fabs(F);
      }
  // Turn force into pressure
  if (pressure)
    for (int i=0; i<size; ++i) std::get<4>(data.at(i)) /= (2*PI*sg[i]);
  // Return data
  return data;
}

inline void Sectorization::firstHalfKick() {
  double dt = 0.5 * epsilon;
  bool doTemp = (temperature>0 && tempDelay<time-lastTemp);
  if (doTemp) {
#pragma vector aligned
#pragma simd
    for (int i=0; i<array_end; ++i) {
      double mass = 1./im[i];
      vx[i] += dt*im[i]*fx[i];
      vy[i] += dt*im[i]*fy[i];
      px[i] += epsilon*vx[i];
      py[i] += epsilon*vy[i];
      fx[i] = gravity.x*mass - 6*PI*viscosity*sg[i]*vx[i]; // Stokes' law
      fy[i] = gravity.y*mass - 6*PI*viscosity*sg[i]*vy[i];
      wrap(px[i], py[i]); // Wrap position
      om[i] += dt*iI[i]*tq[i];
      th[i] += epsilon*om[i];
      wrap(th[i]);        // Wrap theta
      tq[i] = 0;
      // Temperature force
      double DT = DT1*(1./sg[i]); // Assumes Kb = 1;
      double coeffD = sqrt(2*DT)*sqrtTempDelay;
      vec2 TForce = coeffD*randNormal()*randV();
      fx[i] += TForce.x; fy[i] += TForce.y;
    }
  }
  else {
    if (gravity!=Zero) {
      if (!drag) {
#pragma vector aligned
#pragma simd
	for (int i=0; i<array_end; ++i) {
          double mass = 1./im[i];
          vx[i] += dt*im[i]*fx[i];
          vy[i] += dt*im[i]*fy[i];
          px[i] += epsilon*vx[i];
          py[i] += epsilon*vy[i];
          fx[i] = gravity.x*mass;
          fy[i] = gravity.y*mass;
          wrap(px[i], py[i]);
          om[i] += dt*iI[i]*tq[i];
          th[i] += epsilon*om[i];
          tq[i] = 0;
	}
      }
      else {
#pragma vector aligned
#pragma simd
	for (int i=0; i<array_end; ++i) {
	  double mass = 1./im[i];
	  vx[i] += dt*im[i]*fx[i];
	  vy[i] += dt*im[i]*fy[i];
	  px[i] += epsilon*vx[i];
	  py[i] += epsilon*vy[i];
	  fx[i] = gravity.x*mass - 6*PI*viscosity*sg[i]*vx[i]; // Stokes' law
	  fy[i] = gravity.y*mass - 6*PI*viscosity*sg[i]*vy[i];
	  wrap(px[i], py[i]);
	  om[i] += dt*iI[i]*tq[i];
	  th[i] += epsilon*om[i];
	  tq[i] = 0;
	}
      }
    }
    else if (gravity==0) {
      if (!drag) {
#pragma vector aligned
#pragma simd
	for (int i=0; i<array_end; ++i) {
	  double mass = 1./im[i];
	  vx[i] += dt*im[i]*fx[i];
	  vy[i] += dt*im[i]*fy[i];
	  px[i] += epsilon*vx[i];
	  py[i] += epsilon*vy[i];
	  fx[i] = 0;
	  fy[i] = 0;
	  wrap(px[i], py[i]);
	  om[i] += dt*iI[i]*tq[i];
	  th[i] += epsilon*om[i];
	  tq[i] = 0;
	}
      }
      else { // Gravity == 0, drag != 0
#pragma vector aligned
#pragma simd
	for (int i=0; i<array_end; ++i) {
	  double mass = 1./im[i];
	  vx[i] += dt*im[i]*fx[i];
	  vy[i] += dt*im[i]*fy[i];
	  px[i] += epsilon*vx[i];
	  py[i] += epsilon*vy[i];
	  fx[i] = -6*PI*viscosity*sg[i]*vx[i]; // Stokes' law
	  fy[i] = -6*PI*viscosity*sg[i]*vy[i];
	  wrap(px[i], py[i]);
	  om[i] += dt*iI[i]*tq[i];
	  th[i] += epsilon*om[i];
	  tq[i] = 0;
	}
      }
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

inline void Sectorization::wrap(double &x, double &y) {
  if (wrapX) {
    if (x<simBounds.left)       x = simBounds.right-fmod(simBounds.left-x, simBounds.right-simBounds.left);
    else if (simBounds.right<x) x = fmod(x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
  }
  if (wrapY) {
    if (y<simBounds.bottom)     y = simBounds.top-fmod(simBounds.bottom-y, simBounds.top-bounds.bottom);
    else if (simBounds.top<y)   y = fmod(y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
  }
}

inline void Sectorization::wrap(vec2 &pos) {
  wrap(pos.x, pos.y);
}
inline void Sectorization::wrap(double &theta) {
  if (theta<0) theta = 2*PI-fmod(-theta, 2*PI);
  else if (2*PI<theta) theta = fmod(theta, 2*PI);
}

inline int Sectorization::getSec(const vec2 &position) {
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

inline void Sectorization::createNeighborLists(bool force) {
  if (sectors==0) return;
  if (!doInteractions) return;
  
  // Find an upper bound on how far particles might have moved from one another, unless we are forcing a redo
  if (!force) {
    double m0 = 0, m1 = 0;
    for (int i=0; i<array_end; ++i) {
      if (0<=it[i]) {
	double distSqr = sqr(getDisplacement(positionTracker[i], vec2(px[i], py[i])));
	if (distSqr>m0) m0 = distSqr;
	else if (distSqr>m1) m1 = distSqr;
      }
    }
    maxNLDiff = sqrt(m0) + sqrt(m1);
    if (maxNLDiff<skinDepth) return;
  }
  // Update sectors
  updateSectors();
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
	double sigma = sg[i], velocity = sqrt(sqr(vx[i])+sqr(vy[i]));
	double range = particle_cutoff(sigma, it[i]) * max(velocity, 1.);
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
            if (sqr(r)<sqr(range+particle_cutoff(sg[*q],it[*q])+skinDepth)) nlist.push_back(*q);
          }
        // Bottom left
        int sx = x-1, sy = y-1;
	if (0<sy) {
	  if (0<sx) // Don't hit edge sectors
	    for (auto &q : sectors[sy*nsx+sx]) {
	      vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	      if (sqr(r)<sqr(range+particle_cutoff(sg[q],it[i])+skinDepth)) nlist.push_back(q);
	    }
	  // Bottom
	  sx = x;
	  for (auto &q : sectors[sy*nsx+sx]) {
	    vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	    if (sqr(r)<sqr(range+particle_cutoff(sg[q],it[i])+skinDepth)) nlist.push_back(q);
	  }
	}
        // Left
        sx = x-1; sy = y;
	if (0<sx) {
	  for (auto &q : sectors[sy*nsx+sx]) {
	    vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	    if (sqr(r)<sqr(range+particle_cutoff(sg[q],it[i])+skinDepth)) nlist.push_back(q);
	  }
	  // Top left
	  sy = y+1;
	  if (sy<nsy-1)
	    for (auto &q : sectors[sy*nsx+sx]) {
	      vec2 r = getDisplacement(px[i], py[i], px[q], py[q]);
	      if (sqr(r)<sqr(range+particle_cutoff(sg[q],it[i])+skinDepth)) nlist.push_back(q);
	    }
	}
	// Add the neighbor list to the collection if the particle has neighbors
	if (nlist.size()>1) neighborList.push_back(nlist);
      }
    }
  // Update position tracker array
  for (int i=0; i<array_end; ++i) positionTracker[i] = vec2(px[i], py[i]);
}

inline void Sectorization::createEdgeNeighborLists() {
  
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
      if (sqr(displacement)<sqr(1.25*particle_cutoff(sg[i], it[i])))
	lst.push_back(&w);
    }
    if (!lst.empty())
      wallNeighbors.push_back(pair<int, list<Wall*> >(i, lst));
  }
}

inline vec2 Sectorization::getDisplacement(vec2 A, vec2 B) {
  return getDisplacement(A.x, A.y, B.x, B.y);
}

inline vec2 Sectorization::getDisplacement(double ax, double ay, double bx, double by) {
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
  return vec2(X,Y);
}

inline void Sectorization::makeSectors() {
  cutoff = maxCutR+secCutR;
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
}

void Sectorization::updatePList() {
  plist.clear();
  clist.clear();
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
    if (useCharacteristics) clist.push_back(ch[i]);
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
  if (useCharacteristics && ch) delete [] ch;
  // Reallocate
  int tsize = asize + easize;
  pdata[0]  = px = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[1]  = py = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[2]  = vx = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[3]  = vy = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[4]  = fx = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[5]  = fy = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[6]  = th = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[7]  = om = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[8]  = tq = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[9]  = sg = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[10] = im = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[11] = iI = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[12] = rp = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[13] = ds = (double*)aligned_alloc(64, tsize*sizeof(double));
  pdata[14] = cf = (double*)aligned_alloc(64, tsize*sizeof(double));
  it             =    (int*)aligned_alloc(64, tsize*sizeof(int));
  memset(it, -1, tsize*sizeof(int)); // Each int is 4 bytes
  if (useCharacteristics) {
    ch           = (Characteristic**)aligned_alloc(64, asize*sizeof(Characteristic*));
    memset(ch, 0, asize*sizeof(Characteristic*));
  }
  else ch = 0;
  // Set position tracker array
  if (positionTracker) delete [] positionTracker;
  positionTracker = (vec2*)aligned_alloc(64, asize*sizeof(vec2));
  size = 0; esize = 0;
  array_end = 0; earray_end = asize;
}
 
inline void Sectorization::zeroPointers() {
  positionTracker = 0;
  px = py = vx = vy = fx = fy = th = om = tq = sg = im = iI = rp = ds = cf = 0;
  it = 0;
  ch = 0;
  for (int i=0; i<15; ++i) pdata[i] = 0;
  sectors = 0;
}

inline void Sectorization::setParticles() {
  int i=0;
  // Set particles
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
    ++i;
  }
  size = i;
  array_end = i;
  for ( ; i<asize; ++i) it[i] = -1; // No particle stored here
  // Set characteristics
  if (useCharacteristics && ch) {
    i = 0;
    for (auto c : clist) {
      ch[i] = c;
      ++i;
    }
  }
 }

inline void Sectorization::atom_move() {
  // Assumes that only particles that were in boundary might need to move
  //  tl, tm, tr;
  //  ml, **  mr;
  //  bl, bm, br;
  list<int> move_lsts[9]; // Entry 4 will always be unused

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
  // Do the actual migration
  auto start = current_time();
  passParticles(-1, -1, move_lsts[0]); // bl
  passParticles( 0, -1, move_lsts[1]); // bm
  passParticles(+1, -1, move_lsts[2]); // br
  passParticles(-1,  0, move_lsts[3]); // ml
  passParticles(+1,  0, move_lsts[5]); // mr
  passParticles(-1, +1, move_lsts[6]); // tl
  passParticles( 0, +1, move_lsts[7]); // tm
  passParticles(+1, +1, move_lsts[8]); // tr
  auto end = current_time();
  transferTime += time_span(end, start);
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
      buffer[i+15] = static_cast<double>(it[j]);
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
  holes.sort(); // Holes are now in order
  for (auto J=holes.begin(); J!=holes.end(); ++J) {
    int j = *J;
    while (it[array_end-1]<0) array_end--;
    if (array_end-1<=j) break; // No further holes need to be filled
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
    if (ch) {
      ch[j] = ch[array_end-1];
      ch[array_end-1] = 0;
    }
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
  auto start = current_time();
  // Pass left edge to the left
  if (0<dx) // Or wrap y
    for (int y=1; y<nsy-1; ++y)
      for (auto i : sectors[nsx*y+1]) 
	if (0<=it[i]) leftEdge.push_back(i);
  passParticles(-1, 0, leftEdge, true);

  // Pass right edge to the right
  if (dx<ndx-1) // Or wrap x
    for (int y=1; y<nsy-1; ++y)
      for (auto i : sectors[nsx*y+nsx-1]) 
	if (0<=it[i]) leftEdge.push_back(i);
  passParticles(+1, 0, rightEdge, true);

  // Pass top edge upwards
  if (dy<ndy-1)
    for (int x=0; x<nsx; ++x)
      for (auto i : sectors[nsx*(nsy-1)+x]) 
	if (0<=it[i]) topEdge.push_back(i);
  passParticles(0, +1, topEdge, true);

  // Pass bottom edge downwards
  if (0<dy)
    for (int x=0; x<nsx; ++x)
      for(auto i: sectors[nsx+x]) 
	if (0<=it[i]) bottomEdge.push_back(i);
  passParticles(0, -1, bottomEdge, true);
  
  // Create neighborhood lists for the edge sector particles
  createEdgeNeighborLists();
  // Timing
  auto end = current_time();
  transferTime += time_span(end, start);
}

inline void Sectorization::memory_rearrange() {
  if (size==0) return;
  // Get a memory friendly ordering for the particles
  int count = 0;
  vector<int> order(size);
  for (int y=1; y<nsy-1; ++y)
    for (int x=1; x<nsx-1; ++x)
      for (auto index : sectors[nsx*y+x])
	if (-1<it[index]) {
	  order[count] = index;
	  ++count;
	}
  // Create new arrays and put the particles in them in our new order
  double *pd[15]; // New pdata
  int *nit; // New interaction
  // Reallocate
  int tsize = asize + easize;
  pd[0]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[1]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[2]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[3]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[4]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[5]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[6]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[7]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[8]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[9]  = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[10] = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[11] = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[12] = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[13] = (double*)aligned_alloc(64, tsize*sizeof(double));
  pd[14] = (double*)aligned_alloc(64, tsize*sizeof(double));
  nit    =    (int*)   aligned_alloc(64, tsize*sizeof(int));
  // Set
  for (int i=0; i<size; ++i) {
    int j     = order[i];
    pd[0][i]  = px[j];
    pd[1][i]  = py[j];
    pd[2][i]  = vx[j];
    pd[3][i]  = vy[j];
    pd[4][i]  = fx[j];
    pd[5][i]  = fy[j];
    pd[6][i]  = th[j];
    pd[7][i]  = om[j];
    pd[8][i]  = tq[j];
    pd[9][i]  = sg[j];
    pd[10][i] = im[j];
    pd[11][i] = iI[j];
    pd[12][i] = rp[j];
    pd[13][i] = ds[j];
    pd[14][i] = cf[j];
    nit[i]    = it[j];
  }
  // Delete old arrays and set new data
  for (int i=0; i<15; ++i) {
    delete [] pdata[i];
    pdata[i] = pd[i];
  }
  // Set new data arrays
  px = pdata[0]; py = pdata[1]; vx = pdata[2]; vy = pdata[3]; fx = pdata[4]; fy = pdata[5]; th = pdata[6]; om = pdata[7]; tq = pdata[8]; sg = pdata[9]; im = pdata[10]; iI = pdata[11]; rp = pdata[12]; ds = pdata[13]; cf = pdata[14]; it = nit;
  // Rebuild all lists (arrays are already compressed)
  itersSinceBuild = 0;
  createNeighborLists(true); // Position tracker is updated here
  if (doWallNeighbors) createWallNeighborList();
  redoLists = false;
}
