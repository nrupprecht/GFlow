#ifndef __FLOW_HPP__GFLOW__
#define __FLOW_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class Flow : public Modifier {
  public:
    //! @brief Constructor.
    Flow(GFlow*);

    //! @brief Pre forces function, accelerates particles.
    virtual void pre_forces();

  private:
    RealType drag;
  };

}
#endif // __FLOW_HPP__GFLOW__