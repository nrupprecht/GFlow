#include "SimData.hpp"
#include "../data/DataRecord.hpp"
#include "../forces/ExternalForce.hpp"
#include "ForceHandler.hpp"

namespace GFlow {

  SimData::SimData(const Bounds& b, const Bounds& sb) : domain_size(0), domain_capacity(0), edge_size(0), edge_capacity(0), bounds(b), simBounds(sb), wrapX(true), wrapY(true), sectors(nullptr), forceHandler(nullptr) {
    // Set easy access array
    for (int i=0; i<15; ++i) pdata[i] = 0;

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
    // Set summary pointer
    setPData();
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
    walls.push_back(w);
  }

  void SimData::addWall(const Bounds& b) {
    vec2 bl(b.left, b.bottom), br(b.right, b.bottom), tl(b.left, b.top), tr(b.right, b.top);
    walls.push_back(Wall(bl, br));
    walls.push_back(Wall(bl, tl));
    walls.push_back(Wall(tl, tr));
    walls.push_back(Wall(br, tr));
  }

  void SimData::addParticle(const Particle& p) {
    // Check if we have enough room
    if (domain_size==domain_capacity) {
      reserveAdditional(10); // HOW MUCH TO RESERVE?
    }
    // Set arrays
    px[domain_size] = p.position.x;
    py[domain_size] = p.position.y;
    vx[domain_size] = p.velocity.x;
    vy[domain_size] = p.velocity.y;
    fx[domain_size] = p.force.x;
    fy[domain_size] = p.force.y;
    th[domain_size] = p.theta;
    om[domain_size] = p.omega;
    tq[domain_size] = p.torque;
    sg[domain_size] = p.sigma;
    im[domain_size] = p.invMass;
    iI[domain_size] = p.invII;
    rp[domain_size] = p.repulsion;
    ds[domain_size] = p.dissipation;
    cf[domain_size] = p.coeff;
    it[domain_size] = p.interaction;
    // Set position record
    positionRecord[domain_size] = p.position;
    // Increment domain_size
    ++domain_size;
  }

  void SimData::addParticle(const vector<Particle>& particles) {
    for (const auto& p : particles) addParticle(p);
  }

  void SimData::addParticle(const Particle& p, Characteristic *c) {
    characteristics.emplace(domain_size, c);
    addParticle(p);
  }

  RealType** SimData::getPData() {
    return pdata;
  }
  
  vector<Particle> SimData::getParticles() {
    vector<Particle> plist;
    for (int i=0; i<domain_size; ++i) {
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

  void SimData::getPressureData(vector<PData>& positions) {
    // We will sort out particles with it<0 at the end
    vector<PData> pos;
    for (int i=0; i<domain_size; ++i)
      pos.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));

    const auto& verletList = sectors->getVerletList();
    const auto& wallList   = sectors->getWallList();

    // Get data
    forceHandler->pForcesRec(verletList, this, pos);
    forceHandler->wForcesRec(wallList, this, pos);

    // Make sure we only have "real" particles (it>-1)
    for (int i=0; i<domain_size; ++i)
      if (it[i]>-1) {
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
    for (int i=0; i<domain_size; ++i) positionRecord[i] = vec2(px[i], py[i]);
  }

  void SimData::setPData() {
    pdata[0]  = px.getPtr();
    pdata[1]  = py.getPtr();
    pdata[2]  = vx.getPtr();
    pdata[3]  = vy.getPtr();
    pdata[4]  = fx.getPtr();
    pdata[5]  = fy.getPtr();
    pdata[6]  = th.getPtr();
    pdata[7]  = om.getPtr();
    pdata[8]  = tq.getPtr();
    pdata[9]  = sg.getPtr();
    pdata[10] = im.getPtr();
    pdata[11] = iI.getPtr();
    pdata[12] = rp.getPtr();
    pdata[13] = ds.getPtr();
    pdata[14] = cf.getPtr();
  }
  
}
