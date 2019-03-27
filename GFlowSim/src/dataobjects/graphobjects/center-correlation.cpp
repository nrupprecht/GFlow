#include "center-correlation.hpp"
// Other files
#include "../../base/simdata.hpp"
#include "../../utility/vectormath.hpp"

namespace GFlowSimulation {

  CenterCorrelation::CenterCorrelation(GFlow *gflow) : GraphObject(gflow, "CenterCorr", "distance", "counts"), nbins(100), radius(0.15) {
    bins = vector<int>(nbins, 0);
  };

  void CenterCorrelation::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // Get some data
    int size = simData->size();
    RealType dr = radius/nbins;
    Bounds bnds = gflow->getBounds();
    // Create center coordinate
    RealType *center = new RealType[sim_dimensions], *acc = new RealType[sim_dimensions];
    center[0] = 0;
    for (int d=0; d<sim_dimensions; ++d) center[d] = bnds.min[d] + 0.5*bnds.wd(d);
    // Look through all particles
    for (int n=0; n<size; ++n) {
      if (simData->Im(n)>0) { // Only look at movable particles
        // Measure distance from a line at x=0.
        subtractVec(simData->X(n), center, acc, sim_dimensions);
        acc[0] = 0;
        RealType d=dotVec(acc, acc, sim_dimensions);
        // Check if particles are within th radius
        if (d<sqr(radius)) {
          int b = sqrt(d)/dr;
          // Bin data
          if (b<nbins) ++bins[b];
        }
      }
    }
    // Clean up
    delete [] center;
    delete [] acc;
  }

  void CenterCorrelation::setRadius(RealType r) {
    radius = r;
  }

  bool CenterCorrelation::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (bins.empty()) return true;

    // Push data into the data vector
    RealType dr = radius/nbins;
    for (int i=0; i<bins.size(); ++i)
      data.push_back(RPair(i*dr, bins[i]));

    // Graph object does the actual writing.
    return GraphObject::writeToFile(fileName, useName);
  }
}