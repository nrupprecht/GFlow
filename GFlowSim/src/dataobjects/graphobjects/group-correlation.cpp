#include "group-correlation.hpp"
// Other files
#include "../../base/simdata.hpp"
#include "../../utility/vectormath.hpp"

namespace GFlowSimulation {

  GroupCorrelation::GroupCorrelation(GFlow *gflow) : GraphObject(gflow, "GroupCorr", "distance", "counts") {
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
    if (!DataObject::_check()) return;

    // Get arrays
    RealType **x = simData->X();
    int      *id = simData->Id();
    int    *type = simData->Type();
    int     size = simData->size();
    // Find the local ids of the group atoms if simdata has altered the local ids.
    if (locals_changed) {
      local_group.clear();
      for (auto gid : global_group) {
        int local = simData->getLocalID(gid);
        local_group.insert(local);
      }
    }
    locals_changed = false;

    // A lambda for checking if the local id list contains an id.
    auto contains = [&] (int id) -> bool {
      // Try to find the id in the set.
      auto it = local_group.find(id);
      // If the iterator returns an element of local_group, then id is in local group.
      // Therefore, check if it is a valid element, i.e. it != local_group.end().
      return it!=local_group.end();
    };

    // Displacement vector
    RealType *dx = new RealType[sim_dimensions];
    RealType dr  = radius/nbins;

    // Go through all non-group particles
    for (int i=0; i<size; ++i) {
      // Only look at real particles that are not in the group
      if (type[i]<0 || contains(i)) continue;
      // Compute which group particle is closest to particle i.
      RealType minDist = -1; // We can check for this easily
      for (auto id : local_group) {
        // Compute distance
        RealType r = gflow->getDistance(x[id], x[i]);
        // Check if the distance is close enough, and smaller than the min distance so far
        if (r<minDist || minDist==-1) minDist = r;
      }
      // If a close enough group particle was found.
      if (0<=minDist && minDist<=radius) {
        // Compute the bin
        int b = minDist/dr;
        // Bin the data
        ++bins[b];
      }
    }

    // Increment counter
    ++data_iters;

    // Clean up 
    delete [] dx;
  }

  //! \brief Add a particle to the 
  void GroupCorrelation::addToGroup(int g_id) {
    global_group.insert(g_id);
  }

  void GroupCorrelation::setRadius(RealType r) {
    radius = r;
  }

  void GroupCorrelation::setNBins(int nb) {
    nbins = nb;
    bins = vector<int>(nbins, 0);
  }

  int GroupCorrelation::size() {
    return global_group.size();
  }

  bool GroupCorrelation::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (bins.empty() || global_group.empty()) return true;

    // Push data into the data vector
    RealType dr = radius/nbins;
    for (int i=0; i<bins.size(); ++i) {
      RealType r = (i+0.5)*dr;
      data.push_back(RPair(r, 
        static_cast<RealType>(bins[i])/static_cast<RealType>(data_iters))
      );
    }

    // Graph object does the actual writing.
    return GraphObject::writeToFile(fileName, useName);
  }

}