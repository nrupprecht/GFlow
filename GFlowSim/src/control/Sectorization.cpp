#include "Sectorization.hpp"

namespace GFlow {

  Sectorization::Sectorization() : sectors(0), cutoff(0.), skinDepth(0.025), maxCutR(0), secCutR(0), nsx(0), nsy(0), sdx(0), sdy(0), isdx(0), isdy(0), bounds(NullBounds), simBounds(NullBounds), wrapX(true), wrapY(true) {    
    rank = 0;
    numProc = 0;
  }

  Sectorization::~Sectorization() {
    if (sectors) delete [] sectors;
    sectors = 0;
  }

  void Sectorization::initialize(SimData* sd) {
    // Set sim data pointer and bounds
    simData = sd;
    bounds = simData->getBounds();
    simBounds = simData->getSimBounds();
    // Set wrapping
    wrapX = simData->getWrapX();
    wrapY = simData->getWrapY();
    // Find appropriate cutoff radii
    auto plist = simData->getParticles();
    for (const auto &p : plist) {
      double r = p.sigma; // CHANGE THIS TO INCORPORATE INTERACTION TYPE
      if (r>maxCutR) maxCutR = r;
      else if (r>secCutR) secCutR = r;
    }
    // Set up sector array
    makeSectors();
    // Get mpi data
    rank = simData->getRank();
    numProc = simData->getNumProc();
  }

  void Sectorization::sectorize() {
    // Check if there are sectors    
    if (sectors==0) return;
    // Set bounds
    int domain_size = simData->getDomainSize();
    // Get Data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    // Clear out old sectors
    for (int i=0; i<nsx*nsy; ++i) sectors[i].clear();
    // Place in sector 
    for (int i=0; i<domain_size; ++i) {
      int sec_num = getSec(px[i], py[i]);
      sectors[sec_num].push_back(i);
    }
  }

#ifdef USE_MPI

  void Sectorization::atom_move() {
    cout << "Moving atoms.\n";
  }

  void Sectorization::atom_copy() {
    cout << "Copying atoms.\n";
  }

#endif

  void Sectorization::createVerletLists(bool force) {
    // Clear out lists
    verletList.clear();

    // Get position data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    int *it = simData->getItPtr();

    // Check how far the particles have moved 
    if (!force) {
      vec2 *positionRecord = simData->getPRPtr();
      int domain_size = simData->getDomainSize();
      RealType max_moved = 0, sec_moved = 0;
      for (int i=0; i<domain_size; ++i) {
	if (it[i]<0) continue;
	RealType moved = sqr( getDisplacement(px[i], py[i], positionRecord[i].x, positionRecord[i].y) );
	if (moved>max_moved) max_moved = moved;
	else if (moved>sec_moved) sec_moved = moved;
      }
      // Check if a redo may be neccessary
      if (sqrt(max_moved)+sqrt(sec_moved) < skinDepth) return; // No need to redo
    }

    // Create symmetric neighbor lists, look in sectors ( * ) or the sector you are in ( <*> )
    // +---------+
    // | *  x  x |
    // | * <*> x |
    // | *  *  x |
    // +---------+

    // Create new lists
    for (int y=1; y<nsy-1; ++y)
      for (int x=1; x<nsx-1; ++x)
	for (auto p=sec_at(x,y).begin(); p!=sec_at(x,y).end(); ++p) {
	  // Get your data
	  int i = *p;
	  RealType sigma = sg[i];
	  list<int> nlist;
	  nlist.push_back(i); // You are at the head of the list
	  vec2 r;
	  // Sector you are in
	  auto q = p; ++q; // Check only particles ordered after you
	  for (; q!=sec_at(x,y).end(); ++q) {
	    r = getDisplacement(px[i], py[i], px[*q], py[*q]);
	    if (sqr(r) < sqr(sigma + sg[*q] + skinDepth))
	      nlist.push_back(*q);
	  }
	  // Bottom left
	  int sx = x-1, sy = y-1;
	  for (const auto &j : sec_at(sx, sy)) {
	    r = getDisplacement(px[i], py[i], px[j], py[j]);
	    if (sqr(r) < sqr(sigma + sg[j] + skinDepth))
	      nlist.push_back(j);
	  }
	  // Bottom
	  sx = x;
	  for (const auto &j : sec_at(sx, sy)) {
            r = getDisplacement(px[i], py[i], px[j], py[j]);
            if (sqr(r) < sqr(sigma + sg[j] + skinDepth))
              nlist.push_back(j);
          }
	  // Left
	  sx = x-1; sy = y;
	  for (const auto &j : sec_at(sx, sy)) {
            r = getDisplacement(px[i], py[i], px[j], py[j]);
            if (sqr(r) < sqr(sigma + sg[j] + skinDepth))
              nlist.push_back(j);
          }
	  // Top left
	  sy = y+1;
	  for (const auto &j : sec_at(sx, sy)) {
            r = getDisplacement(px[i], py[i], px[j], py[j]);
            if (sqr(r) < sqr(sigma + sg[j] + skinDepth))
              nlist.push_back(j);
          }

	  // Add the neighbor list to the collection if you have neighbors
	  if (nlist.size()>1) verletList.push_back(nlist);
	}
    
    // Update position record
    vec2 *positionRecord = simData->getPRPtr();
    int domain_size = simData->getDomainSize();
    for (int i=0; i<domain_size; ++i) positionRecord[i] = vec2(px[i], py[i]);
    
  }

  void Sectorization::createWallLists(bool force) {
    // Clear out lists
    wallList.clear();

    // Get data
    auto& walls  = simData->getWalls();
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    int *it = simData->getItPtr();
    
    if (!force) {
      vec2 *positionRecord = simData->getPRPtr();
      int domain_size = simData->getDomainSize();
      RealType max_moved = 0, sec_moved = 0;
      for (int i=0; i<domain_size; ++i) {
        if (it[i]<0) continue;
        RealType moved = sqr( getDisplacement(px[i], py[i], positionRecord[i].x, positionRecord[i].y) );
        if (moved>max_moved) max_moved = moved;
      }
      // Check if a redo is neccessary
      if (max_moved < sqr(skinDepth)) return; // No need to redo
    }

    // Create wall list
    int domain_size = simData->getDomainSize();
    for (int i=0; i<walls.size(); ++i) {
      WListSubType lst;
      lst.push_back(i); // Wall is at the head of the list
      for (int j=0; j<domain_size; ++j) {
	if (it[j]<0) continue;
	vec2 displacement = getDisplacement(vec2(px[j], py[j]), walls.at(i).left);
	// Correct the displacement -- THIS DOESNT ALWAYS WORK CORRECTLY
	wallDisplacement(displacement, sg[j], walls.at(i));
	// USE A PARTICLE CUTOFF
	if (sqr(displacement)<sqr(sg[j] + skinDepth)) lst.push_back(j);
      }
      
      if (1<lst.size()) wallList.push_back(lst);
    }
  }

  vec2 Sectorization::getDisplacement(const RealType ax, const RealType ay, const RealType bx, const RealType by) {
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

  vec2 Sectorization::getDisplacement(const vec2 a, const vec2 b) {
    return getDisplacement(a.x, a.y, b.x, b.y);
  }

  inline int Sectorization::getSec(const RealType x, const RealType y) {
    int SX = (x-bounds.left)*isdx+1, SY = (y-bounds.bottom)*isdy+1;
    SX = SX>nsx-1 ? nsx-1 : SX; 
    SX = SX<0 ? 0 : SX;
    SY = SY>nsy-1 ? nsy-1 : SY;
    SY = SY<0 ? 0 : SY;
    return SY*nsx+SX;
  }

  inline void Sectorization::makeSectors() {
    cutoff = maxCutR + secCutR + skinDepth;
    // First estimate
    sdx = sdy = (cutoff+skinDepth);
    nsx = static_cast<int>( max(1., (bounds.right-bounds.left)/sdx) );
    nsy = static_cast<int>( max(1., (bounds.top-bounds.bottom)/sdy) );
    // Actual width and height
    sdx = (bounds.right-bounds.left)/nsx;
    sdy = (bounds.top-bounds.bottom)/nsy;
    isdx = 1./sdx;
    isdy = 1./sdy;
    // Add for edge sectors
    nsx += 2; nsy += 2;
    // Remake sectors
    if (sectors) delete [] sectors;
    sectors = new list<int>[nsx*nsy+1];
  }
  
}
