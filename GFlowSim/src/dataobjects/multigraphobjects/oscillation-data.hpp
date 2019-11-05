#ifndef __OSCILLATION_DATA_HPP__GFLOW__
#define __OSCILLATION_DATA_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"

namespace GFlowSimulation {

  class OscillationData : public MultiGraphObject {
  public:
    //! \brief Constructor
    OscillationData(GFlow*);

    //! \brief Collect the position data from simdata - happens during the post-step phase.
    virtual void post_step() override;

    //! \brief Subtract away average position, so the data is the deviation from the average.
    virtual void post_integrate() override;
  };

}
#endif // __OSCILLATION_DATA_HPP__GFLOW__