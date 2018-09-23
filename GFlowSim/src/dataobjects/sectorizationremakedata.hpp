#ifndef __SECTORIZATION_REMAKE_DATA_HPP__GFLOW__
#define __SECTORIZATION_REMAKE_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  /*
  *  @class SectorizationRemakeData
  *  
  *  Records the number of times sectorization remakes sectors / verlet lists
  *
  */
  class SectorizationRemakeData : public DataObject {
  public:
    // Constructor
    SectorizationRemakeData(GFlow*);

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // The cumulative remakes as a function of time
    vector<RIPair> remakeData;
  };

}
#endif // __SECTORIZATION_REMAKE_DATA_HPP__GFLOW__