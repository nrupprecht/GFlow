#ifndef __POSITION_DATA_HPP__GFLOW__
#define __POSITION_DATA_HPP__GFLOW__

#include "dataobject.hpp"

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
    PositionData(GFlow*, string);
    
    // Destructor
    ~PositionData();

    // Collect the position data from simdata
    virtual void collect();

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    virtual void writeToFile(string, bool=true);

  private:
    // The time steps of when the data was gathered
    vector<RealType> timeStamps;
    // Each time step contains a length ( [number of particles] * [DIMENSIONS] ) array 
    // for the positions of the particles in [DIMENSIONS] dimensions
    vector<RealType*> positions; 
  };

}
#endif // __POSITION_DATA_HPP__GFLOW__