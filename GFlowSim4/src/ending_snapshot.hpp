#ifndef __ENDING_SNAPSHOT_HPP__GFLOW__
#define __ENDING_SNAPSHOT_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  class EndingSnapshot : public DataObject {
  public:
    // Constructor
    EndingSnapshot(GFlow*);
    
    // Destructor
    ~EndingSnapshot();

    // Take a snapshot at the end of the run
    virtual void post_integrate();

    virtual bool writeToFile(string, bool=true) final;

  private:
    // Each time step contains a length ( [number of particles] * [DIMENSIONS] ) array 
    // for the positions of the particles in [DIMENSIONS] dimensions
    RealType* positions; 

    // The amount of data we collect per particle
    int dataWidth;
    
    // The number of particles in the position enties
    int number;
  };

}
#endif // __ENDING_SNAPSHOT_HPP__GFLOW__