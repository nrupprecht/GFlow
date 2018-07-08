#ifndef __MIN_INTERACTING_DISTANCE_HPP__GFLOW__
#define __MIN_INTERACTING_DISTANCE_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  class MinInteractingDistance : public DataObject {
  public:
    // Constructor
    MinInteractingDistance(GFlow*);

    virtual void post_step() final;

    virtual bool writeToFile(string, bool=true) final;

  private:
    vector<RPair> data;
  };

}
#endif // __MIN_INTERACTING_DISTANCE_HPP__GFLOW__