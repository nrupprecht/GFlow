#include "trajectorydata.hpp"

namespace GFlowSimulation {

  TrajectoryData::TrajectoryData(GFlow *gflow) : DataObject(gflow, "Trajectory") {};

  TrajectoryData::~TrajectoryData() {

  }

  void TrajectoryData::post_step() {
    // Data collection:
    
  }

  bool TrajectoryData::writeToFile(string, bool) {
    return false;
  }

}