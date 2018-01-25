#include "../data/DataRecord.hpp"
#include "../forces/ExternalForce.hpp"
#include "ForceHandler.hpp"

template<typename Derived>
SimDataBase<Derived>::SimData(const Bounds& b, const Bounds& sb) : domain_size(0), domain_end(0), domain_capacity(0), edge_size(0), edge_capacity(0), bounds(b), simBounds(sb), wrapX(true), wrapY(true), sectors(nullptr), forceHandler(nullptr), terminate(false) {
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

template<typename Derived>
SimDataBase<Derived>::~SimData() {
  // We are responsible for cleaning up external forces
  for (auto& f : externalForces) delete f;
  // We are responsible for cleaning up characteristics
  for (auto& c : characteristics) delete c.second;
}

template<typename Derived>
void SimDataBase<Derived>::reserve(int dcap, int ecap) {
  this->derived().reserve(dcap, ecap);
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

template<typename Derived>
void SimDataBase<Derived>::reserveAdditional(int d_additional, int e_additional) {
  this->derived().reserveAdditional(d_additional, e_additional);
  /*
  // Check if ecap should be changed
  if (e_additional<0) e_additional = 0;
  // Total room we will need
  int dcap = domain_capacity + d_additional;
  int ecap = edge_capacity + e_additional;
  reserve(dcap, ecap);
  */
}

template<typename Derived>
void SimDataBase<Derived>::addWall(const Wall& w) {
  walls.push_back(w);
}

template<typename Derived>
void SimDataBase<Derived>::addWall(const Bounds& b) {
  vec2 bl(b.left, b.bottom), br(b.right, b.bottom), tl(b.left, b.top), tr(b.right, b.top);
  walls.push_back(Wall(bl, br));
  walls.push_back(Wall(bl, tl));
  walls.push_back(Wall(tl, tr));
  walls.push_back(Wall(br, tr));
}

template<typename Derived>
void SimDataBase<Derived>::addWall(const vector<Wall>& walls) {
  for (const auto &w : walls) addWall(w);
}

template<typename Derived>
int SimDataBase<Derived>::addParticle(const Particle& p) {
  return this->derived().addParticle(p);
  /*
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
  */
}

template<typename Derived>
void SimDataBase<Derived>::addParticle(const vector<Particle>& particles) {
  this->derived().addParticle(particles);
  //for (const auto& p : particles) addParticle(p);
}

template<typename Derived>
int SimDataBase<Derived>::addParticle(const Particle& p, Characteristic *c) {
  this->derived().addParticle(p, c);
  //characteristics.emplace(domain_end, c);
  //return addParticle(p);
}

template<typename Derived>
void SimDataBase<Derived>::removeAt(int index) {
  this->derived().removeAt(i);
  /*
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
  */
}

template<typename Derived>
Particle SimDataBase<Derived>::makeParticle(int i) {
  return this->derived().makeParticle(i);
  /*
    if (domain_end<=i || it[i]<0) throw BadParticle(i);
    // Create and return a copy of the particle
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
    return P;
  */
}

template<typename Derived>
vector<Particle> SimDataBase<Derived>::getParticles() {
  return this->derived().getParticles();
  /*
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
  */
}

template<typename Derived>
vector<Particle> SimDataBase<Derived>::getParticles(vec2 pos, RealType cutoff) {
  return sectors->getParticles(pos, cutoff, this);
}

template<typename Derived>
vector<int> SimDataBase<Derived>::getParticlesID(vec2 pos, RealType cutoff) {
  return sectors->getParticlesID(pos, cutoff, this);
}

template<typename Derived>
void SimDataBase<Derived>::getPressureData(vector<PData>& positions, RealType lowerSizeLimit) {
  /*  
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
  */
}

template<typename Derived>
void SimDataBase<Derived>::wrap(RealType& x, RealType& y) {
  this->derived().wrap(x,y);
  /*
    if (wrapX) {
    if (x<simBounds.left)
    x = simBounds.right-fmod(simBounds.left-x, simBounds.right-simBounds.left);
    else if (simBounds.right<=x)
    x = fmod(x-simBounds.left, simBounds.right-simBounds.left)+simBounds.left;
    }
    if (wrapY) {
    if (y<simBounds.bottom)
    y = simBounds.top-fmod(simBounds.bottom-y, simBounds.top-bounds.bottom);
    else if (simBounds.top<=y)
    y = fmod(y-simBounds.bottom, simBounds.top-simBounds.bottom)+simBounds.bottom;
    }
  */
}

template<typename Derived>
void SimDataBase<Derived>::wrap(RealType& theta) {
  this->derived().wrap(theta);
  //if (theta<0) theta = 2*PI-fmod(-theta, 2*PI);
  //else if (2*PI<=theta) theta = fmod(theta, 2*PI);
}

template<typename Derived>
vec2 SimDataBase<Derived>::getDisplacement(const RealType ax, const RealType ay, const RealType bx, const RealType by) {
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

template<typename Derived>
vec2 SimDataBase<Derived>::getDisplacement(const vec2 a, const vec2 b) {
  return getDisplacement(a.x, a.y, b.x, b.y);
}

template<typename Derived>
vec2 SimDataBase<Derived>::getWallDisplacement(const Wall& w, const vec2 p, RealType sigma) {
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

template<typename Derived>
RealType SimDataBase<Derived>::getPhi() {
  this->derived().getPhi();
  /* 
     RealType vol = 0;
     for (int i=0; i<domain_end; ++i)
     if (-1<it[i]) vol += sqr(sg[i]);
     return PI*vol/simBounds.volume();
  */
}

template<typename Derived>
pair<int, int> SimDataBase<Derived>::getClosestTwo(int id) {
  return sectors->getClosestTwo(id, this);
}

#if USE_MPI == 1
template<typename Derived>
void SimDataBase<Derived>::atomMove() {};

template<typename Derived>
void SimDataBase<Derived>::atomCopy() {};
#endif

template<typename Derived>
void SimDataBase<Derived>::addExternalForce(ExternalForce* force) {
  externalForces.push_back(force);
}

template<typename Derived>
void SimDataBase<Derived>::clearExternalForces() {
  // Clean up the forces
  for (auto& f : externalForces) delete f;
  // Clear the list
  externalForces.clear();
}

template<typename Derived>
void SimDataBase<Derived>::clearValues() {
  this->derived().clearValues();
}

template<typename Derived>
void SimDataBase<Derived>::updatePositionRecord() {
  this->derived().updatePositionRecord();
  //for (int i=0; i<domain_end; ++i) positionRecord[i] = vec2(px[i], py[i]);
}

template<typename Derived>
void SimDataBase<Derived>::setInitialPositions() {
  this->derived().setInitialPositions();
  // initialPositions = vector<vec2>(domain_end);
  // for (int i=0; i<domain_end; ++i) initialPositions[i] = vec2(px[i], py[i]);
}

template<typename Derived>
void SimDataBase<Derived>::removeOverlapping(RealType maxOverlap) {
  if (sectors) sectors->removeOverlapping(maxOverlap);
}

template<typename Derived>
Derived& SimDataBase<Derived>::derived() {
  return static_cast<Derived&>(*this);
}

template<typename Derived>
Derived const& SimDataBase<Derived>::derived() const {
  return static_cast<Derived const&>(*this);
}

template<typename Derived>
inline void SimDataBase<Derived>::compressArrays() {
  // Holes will not in general be in order
  this->derived().compressArrays();
  /*
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
  */
}
