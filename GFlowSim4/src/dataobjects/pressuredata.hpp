#ifndef __PRESSURE_DATA_HPP__GFLOW__
#define __PRESSURE_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class PressureData : public DataObject {
  public:
    //! Constructor
    PressureData(GFlow*);

    //! Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! Write data to a file - if true, the string is a path, and you should use your own name as the file name
    //! Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    //! The vector of data. First entry is time, second entry is the calculated temperature,
    //! the subsequent entries up till the last are the pressures from the different forces. 
    //! The last entry is (P V)/(N T)
    vector<vector<RealType> > data;
  };

}
#endif // __PRESSURE_DATA_HPP__GFLOW__