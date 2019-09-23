#include "group-correlation.hpp"
// Other files
#include "../../base/simdata.hpp"
#include "../../utility/vectormath.hpp"

namespace GFlowSimulation {

  GroupCorrelation::GroupCorrelation(GFlow *gflow) : GraphObject(gflow, "GroupCorr", "distance", "counts"), Group(gflow) {
    bins = vector<int>(nbins, 0);
  };

  void GroupCorrelation::pre_integrate() {
    // Clear data vector
    GraphObject::pre_integrate();
    // Clear bins data
    bins = vector<int>(nbins, 0);
  }

  //! \brief Collect the position data from simdata --- happens during the post-step phase
  void GroupCorrelation::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check() || Group::size()==0) return;

    // Get arrays
    auto x = simData->X();
    auto id = simData->Id();
    auto type = simData->Type();
    auto size = simData->size();
    // Update local ids?
    if (locals_changed) update_local_ids();
    locals_changed = false;

    // A lambda for checking if the local id list contains an id.
    // auto contains = [&] (int id) -> bool {
    //   // Try to find the id in the vector.
    //   auto it = std::find(local_ids.begin(), local_ids.end(), id);
    //   // If the iterator returns an element of local_ids, then id is in local group.
    //   // Therefore, check if it is a valid element, i.e. it != local_ids.end().
    //   return it!=local_ids.end();
    // };

    // Go through all non-group particles
    for (int i=0; i<size; ++i) {
      // Only look at real particles that are not in the group
      if (type[i]<0 || contains(i)) continue;
      // Compute which group particle is closest to particle i.
      RealType minDist = -1; // We can check for this easily
      for (auto id : local_ids) {
        // Compute distance
        RealType r = gflow->getDistance(x[id], x[i]);
        // Check if the distance is close enough, and smaller than the min distance so far
        if (r<minDist || minDist==-1) minDist = r;
      }
      // If a close enough group particle was found.
      if (0<=minDist && minDist<=radius) {
        // Compute the bin
        int b = minDist/bin_width;
        // Bin the data
        ++bins[b];
      }
    }

    // Increment counter
    ++data_iters;
  }

  void GroupCorrelation::setRadius(RealType r) {
    radius = r;
    // Update bin width.
    bin_width = radius/nbins;
  }
  
  void GroupCorrelation::setNBins(int nb) {
    nbins = nb;
    bins = vector<int>(nbins, 0);
    // Update bin width.
    bin_width = radius/nbins;
  }

  bool GroupCorrelation::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (bins.empty() || Group::size()==0) return true;

    // Get the length of the chain.
    RealType group_length = 2*getChainedLength();

    // Push data into the data vector
    for (int i=0; i<bins.size(); ++i) {
      RealType r = (i+0.5)*bin_width;
      data.push_back(RPair(r, 
         static_cast<RealType>(bins[i])/static_cast<RealType>(data_iters) / ((group_length + 2*PI*r) * bin_width)
      ));
    }

    // Graph object does the actual writing.
    return GraphObject::writeToFile(fileName, useName);
  }

}