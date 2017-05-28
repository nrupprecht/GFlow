#include "StandardSectorization.hpp"

namespace GFlow {

  StandardSectorization::StandardSectorization() : Sectorization(), threshold(0) {};

  void StandardSectorization::_createVerletLists() {

    // Get position data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    int *it      = simData->getItPtr();

    // Clear out lists
    verletList.clear();

    // Redo the sectorization so we can make our verlet lists
    _sectorize();

    for (int y=1; y<nsy-1; ++y)
      for (int x=1; x<nsx-1; ++x)
        for (auto p=sec_at(x,y).begin(); p!=sec_at(x,y).end(); ++p) {
          // Get your data
          int i = *p;
          if (it[i]<0) continue;
          RealType sigma = sg[i];
          VListSubType nlist;
          nlist.push_back(i); // You are at the head of the list

	  if (sigma<=threshold) {
	    // Sector you are in
	    auto q = p; ++q; // Check only particles ordered after you
	    for (; q!=sec_at(x,y).end(); ++q) {
	      int j = *q;
	      if (it[j]<0) continue;
	      vec2 r = getDisplacement(px[i], py[i], px[j], py[j]);
	      if (sqr(r) < sqr(sigma + sg[j] + skinDepth))
		nlist.push_back(j);
	    }

	    // Checking lambda
	    auto check = [&] (int sx, int sy) {
	      for (const auto &j : sec_at(sx, sy)) {
		if (it[j]<0) continue;
		if (threshold<sg[j]) continue;
		vec2 r = getDisplacement(px[i], py[i], px[j], py[j]);
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
	  }

	  else { // Scan more than 1 sector around us
	    // Sector sweep
	    RealType scanDist = 2*sigma + skinDepth;
	    int scanBinsX = ceil(scanDist/sdx);
	    int scanBinsY = ceil(scanDist/sdy);
	    int minX = max(1,x-scanBinsX), maxX = min(nsx-2,x+scanBinsX);
	    int minY = max(1,y-scanBinsY), maxY = min(nsy-2,y+scanBinsY);
	    // Sweep
	    for (int sy=minY; sy<=maxY; ++sy)
	      for (int sx=minX; sx<=maxX; ++sx) {
		check(sx,sy,i,it,px,py,sg,sigma,nlist);
	      }
	  }

          // Add the neighbor list to the collection if you have neighbors
          if (nlist.size()>1) verletList.push_back(nlist);
        }

    // Update position record
    simData->updatePositionRecord();
  }

  void StandardSectorization::_makeSectors() {
    // Find a good cutoff
    auto plist = simData->getParticles();
    if (plist.empty()) {
      RealType minSec = 1.;
      nsx = minSec;
      nsy = minSec;
    }
    else {
      // Find the average volume (/ 2 PI) of a particle
      RealType vol = 0;
      for (const auto& p : plist) vol += sqr(p.sigma);
      vol /= plist.size();
      // Find a radius from that volume
      cutoff = 2*sqrt(vol) + skinDepth;
      
      // First estimate of sdx, sdy
      sdx = sdy = cutoff;
      RealType minSec = 1.;
      nsx = static_cast<int>( max(minSec, (bounds.right-bounds.left)/sdx) );
      nsy = static_cast<int>( max(minSec, (bounds.top-bounds.bottom)/sdy) );
    }
    
    // Actual width and height
    sdx = (bounds.right-bounds.left)/nsx;
    sdy = (bounds.top-bounds.bottom)/nsy;
    isdx = 1./sdx;
    isdy = 1./sdy;
    // Threshold radius
    threshold = 0.5*min(sdy-skinDepth, sdy-skinDepth);

    // Add for edge sectors
    nsx += 2; nsy += 2;

    // Remake sectors
    if (sectors) delete [] sectors;
    sectors = new vector<int>[nsx*nsy+1];
  }

  inline void StandardSectorization::check(int sx, int sy, int i, int* it, RealType* px, RealType *py, RealType *sg, RealType sigma, vector<int>& nlist) {
    for (const auto &j : sec_at(sx, sy)) {
      if (it[j]<0) continue;
      if (sg[j]<sigma || px[j]<px[i]) {
	vec2 r = getDisplacement(px[i], py[i], px[j], py[j]);
	if (sqr(r) < sqr(sigma + sg[j] + skinDepth))
	  nlist.push_back(j);
      }
    }
  };

}
