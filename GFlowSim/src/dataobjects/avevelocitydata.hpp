#ifndef __AVE_VELOCITY_DATA_HPP__GFLOW__
#define __AVE_VELOCITY_DATA_HPP__GFLOW__

#include "graphobject.hpp"

namespace GFlowSimulation {

  class AveVelocityData : public GraphObject {
  public:
    //! @brief Constructor.
    AveVelocityData(GFlow*);

    //! @brief Collect the position data from simdata --- happens during the post-step phase.
    virtual void post_step() override;
  };

}
#endif // __AVE_VELOCITY_DATA_HPP__GFLOW__