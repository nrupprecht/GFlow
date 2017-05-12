#include "SimData.hpp"

namespace GFlow {

  SimData::SimData(const Bounds& b, const Bounds& sb) : domain_size(0), domain_capacity(0), edge_size(0), edge_capacity(0), bounds(b), simBounds(sb), wrapX(true), wrapY(true) {};

  void SimData::reserve(int dsize, int esize) {
    // Reserve all data arrays
    px.reserve2(domain_size, dsize, edge_size, esize);
    py.reserve2(domain_size, dsize, edge_size, esize);
    vx.reserve2(domain_size, dsize, edge_size, esize);
    vy.reserve2(domain_size, dsize, edge_size, esize);
    fx.reserve2(domain_size, dsize, edge_size, esize);
    fy.reserve2(domain_size, dsize, edge_size, esize);
    th.reserve2(domain_size, dsize, edge_size, esize);
    om.reserve2(domain_size, dsize, edge_size, esize);
    tq.reserve2(domain_size, dsize, edge_size, esize);
    im.reserve2(domain_size, dsize, edge_size, esize);
    iI.reserve2(domain_size, dsize, edge_size, esize);
    it.reserve2(domain_size, dsize, edge_size, esize);
    // Set sizes
    domain_size = dsize;
    edge_size  = esize;
  }

  void SimData::addWall(const Wall& w) {
    walls.push_back(w);
  }

  void SimData::addParticle(const Particle& p) {
    
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

#ifdef USE_MPI
  void SimData::atomMove() {

  }

  void SimData::atomCopy() {

  }
#endif
  
}
