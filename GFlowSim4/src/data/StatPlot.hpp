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

  typedef void (*StatPlot) (SimData*, vector<RPair>&, const RPair);

  /*
  inline void StatPlot_Velocity(SimData* simData, vector<RPair>& statVector, const RPair bounds) { 
    // Get the bins
    int bins = statVector.size();
    if (bins==0) return;
    // Bin data
    double dq = (bounds.second-bounds.first)/bins;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    int domain_end = simData->getDomainEnd();
    // Gather data
    for (int i=0; i<domain_end; ++i) {
      if (it[i]<0) continue;
      RealType data = sqrt(sqr(vx[i])+sqr(vy[i]));
      int b = (data-bounds.first)/dq;
      if (b<0 || bins<=b) continue; // Keep in bounds
      ++statVector.at(b).second;
    }    
  }

  inline void StatPlot_RadialCorrelation(SimData* simData, vector<RPair>& statVector, const RPair bounds) { 
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
    int domain_end = simData->getDomainEnd();
    int origin = domain_end/2;
    vec2 pos(px[origin], py[origin]);
    // Gather data
    for (int i=0; i<domain_end; ++i) {
      if (it[i]<0) continue;
      RealType data = sqrt(sqr(simData->getDisplacement(pos.x, pos.y, px[i], py[i])));
      if (data==0) continue;
      int b = (data-bounds.first)/dq;
      if (b<0 || bins<=b) continue; // Keep in bounds
      statVector.at(b).second += 1./(2*PI*sqr(data));
    }
  }

  inline void StatPlot_DensityVsDepth(SimData* simData, vector<RPair>& statVector, const RPair bounds) {
    int bins = statVector.size();
    if (bins==0) return;
    // Reset bounds
    RealType bottom = simData->getSimBounds().bottom;
    RealType top   = simData->getSimBounds().top;
    if (statVector.at(0).first==statVector.at(1).first) { // Heights haven't been initialized
      RealType dr = (top-bottom)/bins;
      for (int i=0; i<bins; ++i) statVector.at(i).first = bottom + i*dr;
    }
    // Get data pointers
    RealType *py = simData->getPyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    int domain_end = simData->getDomainEnd();
    RealType dq = (top-bottom)/bins;
    // Bin data
    for (int i=0; i<domain_end; ++i) {
      if (it[i]<0) continue;
      int b = (py[i]-bottom)/dq;
      if (b<0 || bins<=b) continue; // Keep in bounds
      ++statVector.at(b).second;
    }
  }

  inline void StatPlot_PressureVsDepth(SimData* simData, vector<RPair>& statVector, const RPair bounds) {    
    // Get bins
    int bins = statVector.size();
    if (bins==0) return;
    // Reset bounds
    RealType top   = simData->getSimBounds().top;
    RealType bottom = simData->getSimBounds().bottom;
    if (statVector.at(0).first==statVector.at(1).first) { // Heights haven't been initialized
      RealType dr = (top-bottom)/bins;
      for (int i=0; i<bins; ++i) statVector.at(i).first = bottom + i*dr;
    }
    // Gather data
    vector<PData> pos;
    simData->getPressureData(pos);
    // Bin data
    vector<RealType> press(bins, 0);
    vector<int> count(bins, 0);
    double dq = (top-bottom)/bins;
    for (auto& p : pos) {
      int b = (std::get<1>(p)-bottom)/dq;
      if (b<0 || bins<=b) continue;
      press.at(b) += std::get<5>(p);
      ++count.at(b);
    }
    // Combine data
    for (int i=0; i<bins; ++i)
      if (count.at(i)>0) statVector.at(i).second += press.at(i)/count.at(i);
  }

  inline void StatPlot_Alignment(SimData* simData, vector<RPair>& statVector, const RPair bounds) {
    // Get the bins
    int bins = statVector.size();
    if (bins==0) return;
    // Bin data
    double dq = (bounds.second-bounds.first)/bins;
    // Get data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    int domain_end = simData->getDomainEnd();
    // Bin data
    for (int i=0; i<domain_end; ++i) {
      vec2 pos(px[i], py[i]);
      pair<int, int> closest = simData->getClosestTwo(i);
      int j = closest.first, k = closest.second;
      if (j==-1 || k==-1) continue;
      vec2 p1(px[j], py[j]), p2(px[k], py[k]);
      vec2 d1 = simData->getDisplacement(p1,pos), d2 = simData->getDisplacement(p2,pos);
      normalize(d1); normalize(d2);
      RealType theta = acos(d1*d2)/PI; // Cosine of the angle between the closest two particles
      int b = (theta - bounds.first)/dq;
      if (b<0 || bins<=b) ; // Don't go out of bounds
      else statVector.at(b).second += 1;
    }
  }
  */  

}
#endif // __STAT_PLOT_HPP__
