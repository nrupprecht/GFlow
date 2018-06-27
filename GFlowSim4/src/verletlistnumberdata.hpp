#ifndef __VERLET_LIST_NUMBER_DATA_HPP__GFLOW__
#define __VERLET_LIST_NUMBER_DATA_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  /*
  *  @class VerletListNumberData
  *  
  *  Records the number of heads and total particles in verlet lists
  *
  */
  class VerletListNumberData : public DataObject {
  public:
    // Constructor
    VerletListNumberData(GFlow*);

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // The number of heads and total particles in verlet lists
    vector<pair<RealType, IPair> > verletNumbers;
  };

}
#endif // __VERLET_LIST_NUMBER_DATA_HPP__GFLOW__