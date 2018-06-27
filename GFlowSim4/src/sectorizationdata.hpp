#ifndef __SECTORIZATION_DATA_HPP__GFLOW__
#define __SECTORIZATION_DATA_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  class SectorizationData : public DataObject {
  public:
    // Constructor
    SectorizationData(GFlow*);

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // A vector of sector occupations
    vector< vector<int> > sectorOccupation; // [iter] [sector # (linear)] [# of particles]
  };

}
#endif // __SECTORIZATION_DATA_HPP__GFLOW__