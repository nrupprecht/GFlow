#ifndef __VELOCITY_VOLUME_PLOT_HPP__GFLOW__
#define __VELOCITY_VOLUME_PLOT_HPP__GFLOW__

#include "../dataobjecttypes/volumeplotobject2d.hpp"

namespace GFlowSimulation {

  class VelocityVolumePlot : public VolumePlotObject2D {
  public:
    VelocityVolumePlot(GFlow*);

    //! \brief Bin data.
    virtual void post_step() override;
  };

}
#endif // __VELOCITY_VOLUME_PLOT_HPP__GFLOW__