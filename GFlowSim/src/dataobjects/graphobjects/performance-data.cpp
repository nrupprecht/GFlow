#include "performance-data.hpp"

namespace GFlowSimulation {

  PerformanceData::PerformanceData(GFlow *gflow) : GraphObject(gflow, "Performance", "time", "performance") {};

  void PerformanceData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    



  }


}