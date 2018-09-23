#ifndef __AVERAGE_DATA_HPP__GFLOW__
#define __AVERAGE_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class AverageData : public DataObject {
  public:
    // Constructor
    AverageData(GFlow*);

    // Destructor
    ~AverageData();

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    vector<RealType*> data;
    int dataWidth;
  };

}
#endif // __AVERAGE_DATA_HPP__GFLOW__