#ifndef __VERLET_LIST_DATA_HPP__GFLOW__
#define __VERLET_LIST_DATA_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  class VerletListData : public DataObject {
  public:
    // Constructor
    VerletListData(GFlow*);
    
    // Get nForces
    virtual void pre_integrate();

    // Destructor
    ~VerletListData();

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:

    // Records of the verlet lists
    vector<class VerletList> verletLists;

    // Number of forces we recorded verlet lists for
    int nForces;
  };

}
#endif // __VERLET_LIST_DATA_HPP__GFLOW__