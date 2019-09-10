#ifndef __AVERAGE_OMEGA_DATA_HPP__GFLOW__
#define __AVERAGE_OMEGA_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

  class AverageOmegaData : public GraphObject {
  public:
    //! \brief Default constructor.
    AverageOmegaData(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;
  };

}
#endif // __AVERAGE_OMEGA_DATA_HPP__GFLOW__