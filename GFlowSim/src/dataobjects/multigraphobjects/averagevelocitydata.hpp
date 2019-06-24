#ifndef __AVERAGE_VELOCITY_DATA_HPP__GFLOW__
#define __AVERAGE_VELOCITY_DATA_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"

namespace GFlowSimulation {
 
  class AverageVelocityData : public MultiGraphObject {
  public:
    //! @brief Constructor
    AverageVelocityData(GFlow*);

    //! @brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();
  };

}
#endif // __AVERAGE_VELOCITY_DATA_HPP__GFLOW__
