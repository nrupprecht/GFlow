#ifndef __RADIUS_VOLUME_PLOT_HPP__GFLOW__
#define __RADIUS_VOLUME_PLOT_HPP__GFLOW__

#include "../dataobjecttypes/volumeplotobject2d.hpp"

namespace GFlowSimulation {

  class RadiusVolumePlot : public VolumePlotObject2D {
  public:
    //! \brief Constructor
    RadiusVolumePlot(GFlow*);

    //! \brief Bin data.
    virtual void post_step() override;

  private:
    //! \brief We only track particles with radii less than the max_radius.
    RealType max_radius = 1.;

    //! \brief We only track particles with radii greater than the min_radius.
    RealType min_radius = 0.;
  };

}
#endif // __RADIUS_VOLUME_PLOT_HPP__GFLOW__