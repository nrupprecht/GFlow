#ifndef __AVERAGE_SPEED_DATA_HPP__GFLOW__
#define __AVERAGE_SPEED_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class AverageSpeedData : public GraphObject {
  public:
    //! \brief Constructor.
    AverageSpeedData(GFlow*);

    //! \brief Collect the speed data from simdata - happens during the post-step phase.
    virtual void post_step() override;
  };

}
#endif // __AVERAGE_SPEED_DATA_HPP__GFLOW__