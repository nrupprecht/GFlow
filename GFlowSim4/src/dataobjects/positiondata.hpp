#ifndef __POSITION_DATA_HPP__GFLOW__
#define __POSITION_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  /*
  *  @class PositionData
  *  Records the position data of all the objects in the system.
  *  Do this in the post-step phase
  *
  */
  class PositionData : public DataObject {
  public:
    // Constructor
    PositionData(GFlow*);

    // Destructor
    ~PositionData();

    // Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    // Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // The time steps of when the data was gathered
    vector<RealType> timeStamps;

    // Each time step contains a length ( [number of particles] * [DIMENSIONS] ) array 
    // for the positions of the particles in [DIMENSIONS] dimensions
    vector<RealType*> positions; 

    // The amount of data we collect per particle
    int dataWidth;
    
    // The number of particles in the position enties
    vector<int> numbers;
  };

}
#endif // __POSITION_DATA_HPP__GFLOW__