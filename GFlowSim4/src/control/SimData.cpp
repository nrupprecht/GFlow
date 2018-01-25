#include "SimData.hpp"

namespace GFlow {

  SimData::SimData(const Bounds& db, const Bounds& sb) : SimDataBase<SimData>::SimDataBase(db, sb) {};

  const vector<vec2>& SimData::getInitialPositions() const {
    return initialPositions;
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

  int SimData::addParticle(const Particle& p, Characteristic* c) {
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

  Particle SimData::makeParticle(int) {
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

  void SimData::wrap(RealType& x, RealType& y) {
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
  }

  void SimData::wrap(RealType& theta) {
    if (theta<0) theta = 2*PI-fmod(-theta, 2*PI);
    else if (2*PI<=theta) theta = fmod(theta, 2*PI);
  }

  RealType SimData::getPhi() {
    RealType vol = 0;
    for (int i=0; i<domain_end; ++i)
      if (-1<it[i]) vol += sqr(sg[i]);
    return PI*vol/simBounds.volume();
  }

  void SimData::updatePositionRecord() {
    for (int i=0; i<domain_end; ++i) positionRecord[i] = vec2(px[i], py[i]);
  }

  void SimData::setInitialPositions() {
    initialPositions = vector<vec2>(domain_end);
    for (int i=0; i<domain_end; ++i) initialPositions[i] = vec2(px[i], py[i]);
  }

  inline void SimData::compressArrays() {
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
