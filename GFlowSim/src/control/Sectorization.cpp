#include "Sectorization.hpp"
#include "../data/DataRecord.hpp"

namespace GFlow {

  Sectorization::Sectorization() : simData(nullptr), sectors(nullptr), cutoff(0.), skinDepth(default_sectorization_skin_depth), maxCutR(0), secCutR(0), nsx(0), nsy(0), sdx(0), sdy(0), isdx(0), isdy(0), bounds(NullBounds), simBounds(NullBounds), wrapX(true), wrapY(true) {    
    rank = 0;
    numProc = 0;
  }

  Sectorization::Sectorization(SimData* sd) : simData(nullptr), sectors(nullptr), cutoff(0.), skinDepth(default_sectorization_skin_depth), maxCutR(0), secCutR(0), nsx(0), nsy(0), sdx(0), sdy(0), isdx(0), isdy(0), bounds(NullBounds), simBounds(NullBounds), wrapX(true), wrapY(true) {
    rank = 0;
    numProc = 0;
    initialize(sd);
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

    // Set up sector array
    _makeSectors();
    
    // Set this as the sectorization for the sim data
    sd->setSectors(this);

    // Get mpi data
#if USE_MPI == 1
    rank = simData->getRank();
    numProc = simData->getNumProc();
#endif

    // Sectorize
    sectorize();
  }

  void Sectorization::sectorize() {
    // Avoid public virtual functions
    _sectorize();
  }

  void Sectorization::_sectorize() {
    // Check if there are sectors    
    if (sectors==0) return;
    // Set bounds
    int domain_end = simData->getDomainEnd();
    // Get Data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    // Clear out old sectors
    for (int i=0; i<nsx*nsy; ++i) sectors[i].clear();
    // Place in sector 
    for (int i=0; i<domain_end; ++i) {
      int sec_num = getSec(px[i], py[i]);
      sectors[sec_num].push_back(i);
    }
  }

  RealType Sectorization::checkNeedRemake() {
    // Get position data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    int *it = simData->getItPtr();

    // Check how far the particles have moved
    vec2 *positionRecord = simData->getPRPtr();
    int domain_end = simData->getDomainEnd();
    RealType max_moved = 0, sec_moved = 0;
    for (int i=0; i<domain_end; ++i) {
      if (it[i]<0) continue;
      RealType moved = sqr( getDisplacement(px[i], py[i], positionRecord[i].x, positionRecord[i].y) );
      if (moved>max_moved) max_moved = moved;
      else if (moved>sec_moved) sec_moved = moved;
    }
    
    // Calculate max possible movement of particles relative to one another
    movement = sqrt(max_moved)+sqrt(sec_moved);
    // How far have we moved compaired to how far we should move    
    return movement/skinDepth;    
  }

  void Sectorization::createVerletLists() {
    _createVerletLists();
  }

  void Sectorization::_createVerletLists() {
    // Get position data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    int *it      = simData->getItPtr();
    
    // Clear out lists
    verletList.clear();

    // Redo the sectorization so we can make our verlet lists
    _sectorize();

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
	  if (it[i]<0) continue;
	  RealType sigma = sg[i];
	  VListSubType nlist;
	  nlist.push_back(i); // You are at the head of the list
	  vec2 r;
	  // Sector you are in
	  auto q = p; ++q; // Check only particles ordered after you
	  for (; q!=sec_at(x,y).end(); ++q) {
	    int j = *q;
	    if (it[j]<0) continue;
	    r = getDisplacement(px[i], py[i], px[j], py[j]);
	    if (sqr(r) < sqr(sigma + sg[j] + skinDepth)) 
	      nlist.push_back(j);
	  }

	  auto check = [&] (int sx, int sy) {
            for (const auto &j : sec_at(sx, sy)) {
              if (it[j]<0) continue;
              r = getDisplacement(px[i], py[i], px[j], py[j]);
              if (sqr(r) < sqr(sigma + sg[j] + skinDepth))
                nlist.push_back(j);
            }
          };

	  // Bottom left
          int sx = x-1, sy = y-1;
          check(sx, sy);
          // Bottom
          sx = x;
          check(sx, sy);
          // Left
          sx = x-1; sy = y;
          check(sx, sy);
          // Top left
          sy = y+1;
          check(sx, sy);

	  // Add the neighbor list to the collection if you have neighbors
	  if (nlist.size()>1) verletList.push_back(nlist);
	}
    
    // Update position record
    vec2 *positionRecord = simData->getPRPtr();
    int domain_end = simData->getDomainEnd();
    for (int i=0; i<domain_end; ++i) positionRecord[i] = vec2(px[i], py[i]);
  }

  void Sectorization::createWallLists() {
    _createWallLists();
  }

  void Sectorization::_createWallLists() {
    // Get data
    auto& walls  = simData->getWalls();
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    int *it = simData->getItPtr();

    // Only need to do this if there are actually walls
    if (walls.empty()) return;

    // Clear out lists
    wallList.clear();

    // Create wall list
    int domain_end = simData->getDomainEnd();
    for (int i=0; i<walls.size(); ++i) {
      WListSubType lst;
      lst.push_back(i); // Wall is at the head of the list
      for (int j=0; j<domain_end; ++j) {
	if (it[j]<0) continue;
	// vec2 displacement = getDisplacement(vec2(px[j], py[j]), walls.at(i).left);
	// Correct the displacement -- THIS DOESNT ALWAYS WORK CORRECTLY
	// wallDisplacement(displacement, sg[j], walls.at(i));
	vec2 displacement = simData->getWallDisplacement(walls.at(i), vec2(px[j], py[j]), sg[j]);
	// USE A PARTICLE CUTOFF
	if (sqr(displacement)<sqr(sg[j] + skinDepth + default_wall_width)) lst.push_back(j);
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

  int Sectorization::getClosest(int id, SimData *simData) {
    // Get data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    int *it = simData->getItPtr();
    // Search progressively outwards
    vec2 pos(px[id], py[id]);
    int sec = getSec(pos.x, pos.y);
    int sx = sec % nsx, sy = sec/nsx;
    int reach = 1;
    RealType minD = 1e9;
    int minID = -1;
    // Search through sectors for closest particles
    for (int reach=1; reach<min(nsx-2, nsy-2); ++reach) {
      // THIS COULD BE MORE EFFICIENT IF WE ONLY LOOK FOR THINGS AT SQUARE DISTANCE REACH, NOT *WITHIN* REACH
	for (int y=max(sy-reach, 1); y<min(sy+reach, nsy-2); ++y)
	  for (int x=max(sx-reach, 1); x<min(sx+reach, nsx-2); ++x) {
	    RealType threshold = sqr(reach*min(sdx, sdy)); // Corners of the square will be searched before closer regions at the compass points of the square. This keeps us from finding incorrect "closer" particles
	    for (const auto &i : sec_at(x,y)) {
	      if (it[i]<0) continue;
	      vec2 pos2(px[i], py[i]);
	      RealType dsqr = sqr(simData->getDisplacement(pos,pos2));
	      if (dsqr>0 && dsqr<threshold && dsqr<=minD) { // Don't find yourself
		minD = dsqr;
		minID = i;
	      }
	    }
	  }
	// If we have found two, we are done
	if (minID>-1) break;
    }
    return minID;
  }

  pair<int, int> Sectorization::getClosestTwo(int id, SimData *simData) {
    // Get data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    int *it = simData->getItPtr();
    // Search progressively outwards
    vec2 pos(px[id], py[id]);
    int sec = getSec(pos.x, pos.y);
    int sx = sec % nsx, sy = sec / nsx;
    int reach = 1;
    RealType minD1 = 1e9, minD2 = 1e9;
    int min1 = -1, min2 = -1;

    // Helper lambda
    auto check = [&] (int x, int y, RealType threshold) {
      for (const auto &index : sec_at(x,y)) {
	if (it[index]<0) return;
	vec2 pos2(px[index], py[index]);
	RealType dsqr = sqr(simData->getDisplacement(pos,pos2));
	if (dsqr>0 && dsqr<threshold) { // Don't find yourself
	  if (dsqr<=minD1) {
	    minD1 = dsqr;
	    min1 = index;
	  }
	  else if (dsqr<minD2) {
	    minD2 = dsqr;
	    min2 = index;
	  }
	}
      }
    };

    /*
    // Search your sector
    for (const auto &i : sec_at(sx, sy)) check(sx, sy, sqr(min(sdx,sdy)));
    // Start searching surrounding sectors
    for (int reach=1; reach<min(nsx-2, nsy-1); ++reach) {
      RealType threshold = sqr(reach*min(sdx, sdy)); // Corners of the square will be searched before closer regions at the compass points of the square. This keeps us from finding incorrect "closer" particles
      // Left
      if (0<sx-reach)
	for (int y=max(1,sy-reach); y<=min(nsy-2, sy+reach); ++y) check(sx-reach, y, threshold);
      // Top
      if (sy+reach<nsy-1)
	for (int x=max(1,sx-reach); x<=min(nsx-2, sx+reach); ++x) check(x, sy+reach, threshold);
      // Right
      if (sx+reach<nsx-1)
	for (int y=max(1,sy-reach); y<=min(nsy-2, sy+reach); ++y) check(sx+reach, y, threshold);
      // Bottom
      if (0<sy-reach)
	for (int x=max(1,sx-reach); x<=min(nsx-2, sx+reach); ++x) check(x, sy-reach, threshold);
      // If we have found two, we are done
      if (min2>-1) break;
    }
    // Return the pair
    return pair<int, int>(min1, min2);
    */
    
    /*
    // Search through sectors for closest particles
    for (int reach=1; reach<min(nsx-2, nsy-2); ++reach) {
      // THIS COULD BE MORE EFFICIENT IF WE ONLY LOOK FOR THINGS AT SQUARE DISTANCE REACH, NOT *WITHIN* REACH
      for (int y=max(sy-reach, 1); y<min(sy+reach, nsy-2); ++y)
	for (int x=max(sx-reach, 1); x<min(sx+reach, nsx-2); ++x) {
	  RealType threshold = sqr(reach*min(sdx, sdy)); // Corners of the square will be searched before closer regions at the compass points of the square. This keeps us from finding incorrect "closer" particles
	  check(x,y,threshold);
	}
      // If we have found two, we are done
      if (min2>-1) break;
    }
    // Return the pair
    return pair<int, int>(min1, min2);
    */

    // Brute force
    int domain_end = simData->getDomainEnd();
    for (int i=0; i<domain_end; ++i) {
      if (i==id) continue;
      vec2 pos2(px[i], py[i]);
      RealType dsqr = sqr(simData->getDisplacement(pos,pos2));
      if (dsqr < minD1) {
	minD1 = dsqr;
	min1 = i;
      }
      else if (dsqr<minD2) {
	minD2 = dsqr;
	min2 = i;
      }
    }
    return pair<int,int>(min1, min2);
  }

  inline int Sectorization::getSec(const RealType x, const RealType y) {
    int SX = (x-bounds.left)*isdx+1, SY = (y-bounds.bottom)*isdy+1;
    SX = SX>nsx-1 ? nsx-1 : SX; 
    SX = SX<0 ? 0 : SX;
    SY = SY>nsy-1 ? nsy-1 : SY;
    SY = SY<0 ? 0 : SY;
    return SY*nsx+SX;
  }

  void Sectorization::_makeSectors() {
    // Find appropriate cutoff radii
    auto plist = simData->getParticles();
    for (const auto &p : plist) {
      double r = p.sigma; // CHANGE THIS TO INCORPORATE INTERACTION TYPE
      if (r>maxCutR) maxCutR = r;
      else if (r>secCutR) secCutR = r;
    }

    // Cutoff radius
    cutoff = maxCutR + secCutR + skinDepth;

    // First estimate of sdx, sdy
    sdx = sdy = cutoff;
    RealType minSec = 1.;
    nsx = static_cast<int>( max(minSec, (bounds.right-bounds.left)/sdx) );
    nsy = static_cast<int>( max(minSec, (bounds.top-bounds.bottom)/sdy) );

    // Actual width and height
    sdx = (bounds.right-bounds.left)/nsx;
    sdy = (bounds.top-bounds.bottom)/nsy;
    isdx = 1./sdx;
    isdy = 1./sdy;

    // Add for edge sectors
    nsx += 2; nsy += 2;

    // Remake sectors
    if (sectors) delete [] sectors;
    sectors = new vector<int>[nsx*nsy+1];
  }

  // @variable maxOverlap - the overlap threshold
  void Sectorization::removeOverlapping(RealType maxOverlap) {
    // Create verlet lists
    _createVerletLists();
    // Get pointers
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    
    // Set of particles to remove
    std::set<int> remove;

    RealType maxover = 0;

    // Look for overlapping particles
    for (const auto vl : verletList) {
      auto p = vl.begin();
      auto q = p; ++q;
      int i = *p;
      if (remove.find(i)!=remove.end()) continue;
      for (; q!=vl.end(); ++q) {
	int j = *q;
	if (remove.find(j)!=remove.end()) continue;
	vec2 r = simData->getDisplacement(px[i], py[i], px[j], py[j]);
	RealType overlap = sg[i] + sg[j] - sqrt(sqr(r)); // > 0 -> overlap
	RealType over = max(overlap/sg[i], overlap/sg[j]);
	// Mark for removal if neccessary
	if (over>maxOverlap) {
	  if (sg[i]>sg[j]) remove.emplace(j);
	  else             remove.emplace(i);
	}
      }
    }
    // Found all the particles to remove
    for (auto p : remove) simData->removeAt(p);
  }
  
  vector<Particle> Sectorization::getParticles(vec2 pos, RealType scanDist, SimData *simData) {
    // Sectorize
    _sectorize();
    // List of particles to fill
    vector<Particle> particles;
    // Set up
    int sec = getSec(pos.x, pos.y);
    int x = sec%nsx, y = sec/nsx;
    RealType scanDistSqr = sqr(scanDist);
    int scanBinsX = ceil(scanDist/sdx);
    int scanBinsY = ceil(scanDist/sdy);
    int minX = x-scanBinsX, maxX = x+scanBinsX;
    int minY = y-scanBinsY, maxY = y+scanBinsY;
    // Adjust for wrapping or bounds
    if (minX<1) minX = (wrapX ? minX+nsx-2 : 1);
    if (nsx-2<maxX) maxX = (wrapX ? maxX-nsx+2 : nsx-2);
    if (minY<1) minY = (wrapY ? minY+nsy-2 : 1);
    if (nsy-2<maxY) maxY = (wrapY ? maxY-nsy+2 : nsy-2);
    // Sweep
    RealType *px = simData->getPxPtr(), *py = simData->getPyPtr();
    int *it = simData->getItPtr();
    for (int sy = minY; sy!=maxY+1; ++sy) {
      for (int sx = minX; sx!=maxX+1; ++sx) {
	for (const auto I : sec_at(sx, sy)) {
	  if (-1<it[I] && sqr(simData->getDisplacement(pos.x, pos.y, px[I], py[I]))<scanDistSqr)
	    particles.push_back( simData->makeParticle(I) );
	}
	if (wrapX && nsx-1==sx) sx = 0; // The ++ will make this 1
      }
      if (wrapY && nsy-1==sy) sy = 0; // The ++ will make this 1
    }
    return particles;
  }

}
