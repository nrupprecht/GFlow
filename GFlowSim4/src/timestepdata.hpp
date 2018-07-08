#ifndef __TIME_STEP_DATA_HPP__GFLOW__
#define __TIME_STEP_DATA_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  class TimeStepData : public DataObject {
  public:
    // Constructor
    TimeStepData(GFlow*);

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // The data
    vector<RPair> data;
  };

}
#endif // __TIME_STEP_DATA_HPP__GFLOW__