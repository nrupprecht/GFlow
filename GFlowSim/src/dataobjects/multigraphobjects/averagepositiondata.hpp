#ifndef __AVERAGE_POSITION_DATA_HPP__GFLOW__
#define __AVERAGE_POSITION_DATA_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"

namespace GFlowSimulation {
 
  class AveragePositionData : public MultiGraphObject {
  public:
    //! @brief Constructor
    AveragePositionData(GFlow*);

    //! @brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

  private:
    //! \brief Whether to record average position, or sum of positions.
    bool useAve = true;
  };

}
#endif // __AVERAGE_POSITION_DATA_HPP__GFLOW__
