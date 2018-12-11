#ifndef __AVERAGE_POSITION_DATA_HPP__GFLOW__
#define __AVERAGE_POSITION_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {
 
  class AveragePositionData : public DataObject {
  public:
    //! @brief Constructor
    AveragePositionData(GFlow*);

    //! @brief Destructor
    ~AveragePositionData();

    //! @brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! @brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    //! @brief The data
    vector<RPair> *posdata;
  };

}
#endif // __AVERAGE_POSITION_DATA_HPP__GFLOW__
