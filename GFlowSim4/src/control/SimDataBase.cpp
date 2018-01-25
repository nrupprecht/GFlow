#include "../data/DataRecord.hpp"
#include "../forces/ExternalForce.hpp"
#include "ForceHandler.hpp"

template<typename Derived>
SimDataBase<Derived>::SimDataBase(const Bounds& b, const Bounds& sb) : domain_size(0), domain_end(0), domain_capacity(0), edge_size(0), edge_capacity(0), bounds(b), simBounds(sb), wrapX(true), wrapY(true), sectors(nullptr), forceHandler(nullptr), terminate(false) {
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
SimDataBase<Derived>::~SimDataBase() {
  // We are responsible for cleaning up external forces
  for (auto& f : externalForces) delete f;
  // We are responsible for cleaning up characteristics
  for (auto& c : characteristics) delete c.second;
}

template<typename Derived>
void SimDataBase<Derived>::reserve(int dcap, int ecap) {
  this->derived().reserve(dcap, ecap);
}

template<typename Derived>
void SimDataBase<Derived>::reserveAdditional(int d_additional, int e_additional) {
  this->derived().reserveAdditional(d_additional, e_additional);
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
}

template<typename Derived>
void SimDataBase<Derived>::addParticle(const vector<Particle>& particles) {
  this->derived().addParticle(particles);
}

template<typename Derived>
int SimDataBase<Derived>::addParticle(const Particle& p, Characteristic *c) {
  this->derived().addParticle(p, c);
}

template<typename Derived>
void SimDataBase<Derived>::removeAt(int index) {
  this->derived().removeAt(i);
}

template<typename Derived>
Particle SimDataBase<Derived>::makeParticle(int i) {
  return this->derived().makeParticle(i);
}

template<typename Derived>
vector<Particle> SimDataBase<Derived>::getParticles() {
  return this->derived().getParticles();
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
}

template<typename Derived>
void SimDataBase<Derived>::wrap(RealType& theta) {
  this->derived().wrap(theta);
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
}

template<typename Derived>
void SimDataBase<Derived>::setInitialPositions() {
  this->derived().setInitialPositions();
}

template<typename Derived>
const vector<vec2>& SimDataBase<Derived>::getInitialPositions() const {
  return this->derived().getInitialPositions();
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
}
