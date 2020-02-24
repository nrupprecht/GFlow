#ifndef __NUMBER_DATA_HPP__GFLOW__
#define __NUMBER_DATA_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"

namespace GFlowSimulation {

  class NumberData : public MultiGraphObject {
  public:
    //! \brief Default constructor.
    NumberData(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase.
    virtual void post_step() override;
  };

}
#endif // __NUMBER_DATA_HPP__GFLOW__