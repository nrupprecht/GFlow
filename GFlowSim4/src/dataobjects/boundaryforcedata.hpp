#ifndef __BOUNDARY_FORCE_DATA_HPP__GFLOW__
#define __BOUNDARY_FORCE_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class BoundaryForceData : public DataObject {
  public:
    // Constuctor 
    BoundaryForceData(GFlow*);

    void post_step();

    bool writeToFile(string, bool);

    RealType getAverage() const;
  private:
    vector<RPair> bForces;
  };

}
#endif // __BOUNDARY_FORCE_DATA_HPP__GFLOW__