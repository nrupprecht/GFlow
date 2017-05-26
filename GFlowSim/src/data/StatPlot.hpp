/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 25, 2017
 *
 */

#ifndef __STAT_PLOT_HPP__
#define __STAT_PLOT_HPP__

// Includes
#include "../control/SimData.hpp"

namespace GFlow {

  typedef void (*StatPlot) (SimData*, vector<vec2>&, const RPair);

  inline void StatPlot_Velocity(SimData* simData, vector<vec2>& statVector, const RPair bounds) { 
    int bins = statVector.size();
    if (bins==0) return;
    // Bin data
    double dq = (bounds.second-bounds.first)/bins;

    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    int domain_size = simData->getDomainSize();
    // Gather data
    for (int i=0; i<domain_size; ++i) {
      if (it[i]<0) continue;
      RealType data = sqrt(sqr(vx[i])+sqr(vy[i]));
      int b = (data-bounds.first)/dq;
      if (b<0 || bins<=b) continue; // Keep in bounds
      ++statVector.at(b).y;
    }    
  }

  inline void StatPlot_RadialCorrelation(SimData* simData, vector<vec2>& statVector, const RPair bounds) { 
    int bins = statVector.size();
    if (bins==0) return;
    // Bin data
    RealType dq = (bounds.second-bounds.first)/bins;
    // Get data pointers
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Choose an origin particle
    int domain_size = simData->getDomainSize();
    int origin = domain_size/2;
    vec2 pos(px[origin], py[origin]);
    // Gather data
    for (int i=0; i<domain_size; ++i) {
      if (it[i]<0) continue;
      RealType data = sqrt(sqr(simData->getDisplacement(pos.x, pos.y, px[i], py[i])));
      if (data==0) continue;
      int b = (data-bounds.first)/dq;
      if (b<0 || bins<=b) continue; // Keep in bounds
      statVector.at(b).y += 1./(2*PI*sqr(data));
    }
  }

  inline void StatPlot_DensityVsDepth(SimData* simData, vector<vec2>& statVector, const RPair bounds) {
    int bins = statVector.size();
    if (bins==0) return;
    // Reset bounds
    RealType bottom = simData->getSimBounds().bottom;
    RealType top   = simData->getSimBounds().top;
    if (statVector.at(0).x==statVector.at(1).x) { // Heights haven't been initialized
      RealType dr = (top-bottom)/bins;
      for (int i=0; i<bins; ++i) statVector.at(i).x = bottom + i*dr;
    }
    // Get data pointers
    RealType *py = simData->getPyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    int domain_size = simData->getDomainSize();
    RealType dq = (top-bottom)/bins;
    // Bin data
    for (int i=0; i<domain_size; ++i) {
      if (it[i]<0) continue;
      int b = (py[i]-bottom)/dq;
      if (b<0 || bins<=b) continue; // Keep in bounds
      ++statVector.at(b).y;
    }
  }
  
}
#endif // __STAT_PLOT_HPP__
