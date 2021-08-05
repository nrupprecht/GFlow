#include <dataobjects/graphobjects/center-correlation.hpp>
// Other files
#include <base/simdata.hpp>
#include <utility/vectormath.hpp>

using namespace GFlowSimulation;

CenterCorrelation::CenterCorrelation(GFlow *gflow)
    : GraphObject(gflow, "CenterCorr", "distance", "counts"), radius(0.15), nbins(100) {
  bins = vector<int>(nbins, 0);
};

void CenterCorrelation::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }
  // Get some data
  RealType dr = radius / nbins;
  Bounds bnds = gflow->getBounds();
  // Create center coordinate
  Vec center(sim_dimensions), acc(sim_dimensions);
  center[0] = 0;
  for (int d = 0; d < sim_dimensions; ++d) {
    center[d] = bnds.min[d] + 0.5 * bnds.wd(d);
  }
  // Look through all particles
  for (int n = 0; n < simData->size_owned(); ++n) {
    if (simData->Im(n) > 0) { // Only look at movable particles
      // Measure distance from a line at x=0.
      subtractVec(simData->X(n), center.data, acc.data, sim_dimensions);
      acc[0] = 0;
      RealType d = acc * acc;
      // Check if particles are within th radius
      if (d < sqr(radius)) {
        int b = std::sqrt(d) / dr;
        // Bin data
        if (b < nbins) {
          ++bins[b];
        }
      }
    }
  }
}

void CenterCorrelation::setRadius(RealType r) {
  radius = r;
}

bool CenterCorrelation::writeToFile(string fileName, bool useName) {
  // Check if there's anything to do
  if (bins.empty()) {
    return true;
  }

  // Push data into the data vector
  RealType dr = radius / nbins;
  for (int i = 0; i < bins.size(); ++i) {
    data.emplace_back(i * dr, bins[i]);
  }

  // Graph object does the actual writing.
  return GraphObject::writeToFile(fileName, useName);
}
