#include "SimData.hpp"
#include "../data/DataRecord.hpp"
#include "../forces/ExternalForce.hpp"
#include "ForceHandler.hpp"

namespace GFlow {

  SimData::SimData(const Bounds& b, const Bounds& sb) : domain_size(0), domain_end(0), domain_capacity(0), edge_size(0), edge_capacity(0), bounds(b), simBounds(sb), wrapX(true), wrapY(true), sectors(nullptr), forceHandler(nullptr), terminate(false) {
    // Set up MPI (possibly)
#if USE_MPI == 1
#if _CLANG_ == 1
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
#else
    rank = MPI::COMM_WORLD.Get_rank();
    numProc = MPI::COMM_WORLD.Get_size();
#endif
#endif
  };
  
  SimData::~SimData() {
    // We are responsible for cleaning up external forces
    for (auto& f : externalForces) delete f;
    // We are responsible for cleaning up characteristics
    for (auto& c : characteristics) delete c.second;
  }

  void SimData::reserve(int dcap, int ecap) {
    // Check if ecap should be changed
    if (ecap<0) ecap = edge_capacity;
    // Reserve all data arrays
    px.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    py.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    vx.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    vy.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    fx.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    fy.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    th.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    om.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    tq.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    sg.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    im.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    iI.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    rp.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    ds.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    cf.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    it.reserve2(domain_capacity, dcap, edge_capacity, ecap, -1);
    // Reserve positions record
    positionRecord.reserve2(domain_capacity, dcap, edge_capacity, ecap);
    // Set sizes
    domain_capacity = dcap;
    edge_capacity  = ecap;
  }

  void SimData::reserveAdditional(int d_additional, int e_additional) {
    // Check if ecap should be changed
    if (e_additional<0) e_additional = 0;
    // Total room we will need
    int dcap = domain_capacity + d_additional;
    int ecap = edge_capacity + e_additional;
    reserve(dcap, ecap);
  }

  void SimData::addWall(const Wall& w) {
    vec2 left = w.left, right = w.getRight();
    RealType EPS = 0.001;
    // Clip walls so they are always within the proper region if we are wrapping
    if (wrapX || wrapY) {
      if (left.x==right.x) { // Verticle wall - cannot parameterize by x
	// Only need to add (at most) one wall
	if (wrapY) {
	  if (right.y<simBounds.bottom || simBounds.top<left.y) return; // Wall is out of bounds
	  else {
	    left.y = max(left.y, simBounds.bottom+EPS);
	    right.y = min(right.y, simBounds.top-EPS);
	  }
	}
	// Make sure we are not on an edge if we are wrapping
	if (wrapX && left.x==simBounds.left) {
	  left.x += EPS; right.x += EPS;
	}
	if (wrapX && left.x==simBounds.right) {
	  left.x -= EPS; right.x -= EPS;
	}
	// Add the wall
	walls.push_back(Wall(left, right)); // No need to clip
      }
      else { // Parameterize by x
	// SIMPLE CLIP
	// Wall slope
	RealType slope = (right.y-left.y)/(right.x-left.x);
	// Find first x value where the wall is right of the left wall
      }
    }
    else walls.push_back(w);
  }

  void SimData::addWall(const Bounds& b) {
    vec2 bl(b.left, b.bottom), br(b.right, b.bottom), tl(b.left, b.top), tr(b.right, b.top);
    walls.push_back(Wall(bl, br));
    walls.push_back(Wall(bl, tl));
    walls.push_back(Wall(tl, tr));
    walls.push_back(Wall(br, tr));
  }

  int SimData::addParticle(const Particle& p) {
    // Check if we have enough room
    if (domain_end==domain_capacity) {
      reserveAdditional(10); // HOW MUCH TO RESERVE?
    }
    // Set arrays
    px[domain_end] = p.position.x;
    py[domain_end] = p.position.y;
    vx[domain_end] = p.velocity.x;
    vy[domain_end] = p.velocity.y;
    fx[domain_end] = p.force.x;
    fy[domain_end] = p.force.y;
    th[domain_end] = p.theta;
    om[domain_end] = p.omega;
    tq[domain_end] = p.torque;
    sg[domain_end] = p.sigma;
    im[domain_end] = p.invMass;
    iI[domain_end] = p.invII;
    rp[domain_end] = p.repulsion;
    ds[domain_end] = p.dissipation;
    cf[domain_end] = p.coeff;
    it[domain_end] = p.interaction;
    // Set position record
    positionRecord[domain_end] = p.position;
    // Increment domain_size
    ++domain_size;
    ++domain_end;
    // Return where we placed the particle
    return domain_end-1;
  }

  void SimData::addParticle(const vector<Particle>& particles) {
    for (const auto& p : particles) addParticle(p);
  }

  int SimData::addParticle(const Particle& p, Characteristic *c) {
    characteristics.emplace(domain_end, c);
    return addParticle(p);
  }

  void SimData::removeAt(int index) {
    if (it[index]<0) return;
    it[index] = -1;
    auto c = characteristics.find(index);
    if (c!=characteristics.end()) {
      delete c->second;
      characteristics.erase(index);
    }
    --domain_size;
    if (index==domain_end-1) --domain_end;
    holes.push_back(index);
  }
  
  vector<Particle> SimData::getParticles() {
    vector<Particle> plist;
    for (int i=0; i<domain_end; ++i) {
      if (it.at(i)>-1) {
	Particle P(px.at(i), py.at(i), sg.at(i));
	P.velocity    = vec2(vx.at(i), vy.at(i));
	P.force       = vec2(fx.at(i), fy.at(i));
	P.theta       = th.at(i);
	P.omega       = om.at(i);
	P.torque      = tq.at(i);
	P.invMass     = im.at(i);
	P.invII       = iI.at(i);
	P.repulsion   = rp.at(i);
	P.dissipation = ds.at(i);
	P.coeff       = cf.at(i);
	P.interaction = it.at(i);
	// Push particle into list
	plist.push_back(P);
      }
    }
    return plist;
  }

  void SimData::getPressureData(vector<PData>& positions, RealType lowerSizeLimit) {
    // We will sort out particles with it<0 at the end
    vector<PData> pos;
    for (int i=0; i<domain_end; ++i)
      pos.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));

    const auto& verletList = sectors->getVerletList();
    const auto& wallList   = sectors->getWallList();

    // Get data
    forceHandler->pForcesRec(verletList, this, pos);
    forceHandler->wForcesRec(wallList, this, pos);

    // Make sure we only have "real" particles (it>-1)
    for (int i=0; i<domain_end; ++i)
      if (it[i]>-1 && lowerSizeLimit<sg[i]) {
        PData pdata = pos.at(i);
	std::get<5>(pdata) /= 2*PI*sg[i]; // Convert to pressure
        positions.push_back(pdata);
      }
  }

  void SimData::wrap(RealType& x, RealType& y) {
    if (wrapX) {
      if (x<simBounds.left)
	x = simBounds.right-fmod(simBounds.left-x, simBounds.right-simBounds.left);
      else if (simBounds.right<x) 
	x = fmod(x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
    }
    if (wrapY) {
      if (y<simBounds.bottom)
	y = simBounds.top-fmod(simBounds.bottom-y, simBounds.top-bounds.bottom);
      else if (simBounds.top<y)
	y = fmod(y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
    }
  }

  void SimData::wrap(RealType& theta) {
    if (theta<0) theta = 2*PI-fmod(-theta, 2*PI);
    else if (2*PI<theta) theta = fmod(theta, 2*PI);
  }

  vec2 SimData::getDisplacement(const RealType ax, const RealType ay, const RealType bx, const RealType by) {
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

  vec2 SimData::getDisplacement(const vec2 a, const vec2 b) {
    return getDisplacement(a.x, a.y, b.x, b.y);
  }

  vec2 SimData::getWallDisplacement(const Wall& w, const vec2 p, RealType sigma) {
    if (wrapX || wrapY) {
      // Displacement
      vec2 d1 = p-w.left;
      wallDisplacement(d1, sigma, w);
      RealType dsqr1 = sqr(d1);
      // Wrapped displacement
      vec2 d2 = getDisplacement(p, w.left);
      wallDisplacement(d2, sigma, w);
      RealType dsqr2 = sqr(d2);
      // Which is smaller?
      return dsqr1<dsqr2 ? d1 : d2;
    }
    else {
      vec2 d = p-w.left;
      wallDisplacement(d, sigma, w);
      return d;
    }
  }

  RealType SimData::getPhi() {
    RealType vol = 0;
    for (int i=0; i<domain_end; ++i)
      if (-1<it[i]) vol += sqr(sg[i]);
    return PI*vol/simBounds.volume();
  }

#if USE_MPI == 1
  void SimData::atomMove() {
    
  }

  void SimData::atomCopy() {

  }
#endif

  void SimData::addExternalForce(ExternalForce* force) {
    externalForces.push_back(force);
  }

  void SimData::clearExternalForces() {
    // Clean up the forces
    for (auto& f : externalForces) delete f;
    // Clear the list
    externalForces.clear();
  }

  void SimData::clearForceTorque() {
    fx.clearValues();
    fy.clearValues();
    tq.clearValues();
  }

  void SimData::updatePositionRecord() {
    for (int i=0; i<domain_end; ++i) positionRecord[i] = vec2(px[i], py[i]);
  }

  void SimData::removeOverlapping(RealType maxOverlap) {
    if (sectors) sectors->removeOverlapping(maxOverlap);
  }

  inline void SimData::compressArrays() {
    // Holes will not in general be in order
    if (holes.empty()) return;
    holes.sort();
    for (auto J=holes.begin(); J!=holes.end(); ++J) {
      int j = *J;
      while (it[domain_end-1]<0) domain_end--;
      if (domain_end-1<=j) break; // No further holes need to be filled
      px[j] = px[domain_end-1];
      py[j] = py[domain_end-1];
      vx[j] = vx[domain_end-1];
      vy[j] = vy[domain_end-1];
      fx[j] = fx[domain_end-1];
      fy[j] = fy[domain_end-1];
      th[j] = th[domain_end-1];
      om[j] = om[domain_end-1];
      tq[j] = tq[domain_end-1];
      sg[j] = sg[domain_end-1];
      im[j] = im[domain_end-1];
      iI[j] = iI[domain_end-1];
      rp[j] = rp[domain_end-1];
      ds[j] = ds[domain_end-1];
      cf[j] = cf[domain_end-1];
      it[j] = it[domain_end-1];
      
      auto c = characteristics.find(domain_end-1);
      if (c!=characteristics.end()) characteristics.emplace(j, c->second);
      // This entry is now empty
      it[domain_end-1] = -1;    
      // Increment array end
      --domain_end;
    }
    holes.clear();
  }
  
}
