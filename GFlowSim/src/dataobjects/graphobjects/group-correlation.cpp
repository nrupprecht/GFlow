#include "group-correlation.hpp"
// Other files
#include "../../base/simdata.hpp"
#include "../../utility/vectormath.hpp"
#include "../../base/interactionhandler.hpp"

namespace GFlowSimulation {

  GroupCorrelation::GroupCorrelation(GFlow *gflow) : GraphObject(gflow, "GroupCorr", "distance", "counts"), Group(gflow) {
    bins = vector<int>(nbins, 0);
  };

  void GroupCorrelation::pre_integrate() {
    // Clear data vector
    GraphObject::pre_integrate();
    // Clear bins data
    bins = vector<int>(nbins, 0);
    average_length = 0;
    data_iters = 0;
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
    
    RealType minDist = -1; // We can check for this easily
    const RealType margin = 1.5;

    // For finding distance to chain.
    Vec normal(sim_dimensions), dx(sim_dimensions);
    // A function for correcting distances.
    auto check_adjacent = [&] (int i, int id_min, int id2) {
      // Get the vectors we need.
      gflow->getDisplacement(x(i), x(id_min), dx.data);
      gflow->getDisplacement(x(id2), x(id_min), normal.data);
      normal.normalize();
      // Project the displacement from the closest particle to particle i onto the chain from the closest particle
      // to its neighbor in the chain.
      RealType project = normal*dx;
      // If the outside particle is closer to the chain than to the group particle, correct the distance.
      if (project>0) {
        dx -= project*normal;
        RealType new_distance = magnitude(dx);
        if (new_distance<minDist) minDist = new_distance;
      }
    };

    // Look for particles near the chain using the interaction handler's get all within function.    
    vector<int> neighbors;
    std::set<int> all_neighbors;
    for (int i=0; i<Group::size(); ++i) {
      handler->getAllWithin(Group::at(i), neighbors, margin*radius);
      for (auto id : neighbors) all_neighbors.insert(id);
    }

    // Check which of the potential neighbors are close enough.
    for (auto i : all_neighbors) {
      // Only look at real particles that are not in the group
      if (type[i]<0 || contains(i)) continue;
      // Compute which group particle is closest to particle i.
      minDist = -1; // We can check for this easily
      int index = 0, min_index = -1; // Keep track of the index (in the group) of the close particle.
      for (auto id : local_ids) {
        // Compute distance
        RealType r = gflow->getDistance(x(id), x(i));
        // Check if the distance is close enough, and smaller than the min distance so far
        if (r<minDist || minDist==-1) {
          minDist = r;
          min_index = index;
        }
        ++index;
      }      
      // Correct if we want the distance-from-chain.
      if (find_chain_distance && minDist<margin*radius) {
        // Check left segement, if one exists: X---(min)
        if (min_index>0) check_adjacent(i, local_ids[min_index], local_ids[min_index-1]);
        // Check right segment, if one exists: (min)---X 
        if (min_index+1<Group::size()) check_adjacent(i, local_ids[min_index], local_ids[min_index+1]);
      }

      // If a close enough group particle was found.
      if (minDist<=radius) {
        // Compute the bin
        int b = minDist/bin_width;
        // Bin the data
        ++bins[b];
      }
    }

    // Update average length.
    average_length += getChainedLength();
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
    average_length /= data_iters;
    RealType group_length = 2*average_length;
    RealType rho = (simData->number() - Group::size()) / gflow->getBounds().vol(); // Number density

    // Push data into the data vector
    for (int i=0; i<bins.size(); ++i) {
      RealType r = (i+0.5)*bin_width;
      data.push_back(RPair(r, 
         static_cast<RealType>(bins[i])/static_cast<RealType>(data_iters) / ((group_length + 2*PI*r) * bin_width * rho)
      ));
    }

    // Graph object does the actual writing.
    return GraphObject::writeToFile(fileName, useName);
  }

}