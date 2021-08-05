#ifndef __PHI_DATA_HPP__GFLOW__
#define __PHI_DATA_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"

namespace GFlowSimulation {

  class PhiData : public MultiGraphObject {
  public:
    //! \brief Default constructor.
    PhiData(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase.
    virtual void post_step() override;

  private:

    //! \brief The total amount of volume that should be excluded whan calculating the total volume
    //! of the simulation.
    real total_excluded_volume = 0.f;
  };

}
#endif // __PHI_DATA_HPP__GFLOW__