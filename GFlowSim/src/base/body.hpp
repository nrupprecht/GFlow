#ifndef __BODY_HPP__GFLOW__
#define __BODY_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class Body : public Base {
  public:
    //! \brief Default constructor.
    Body(GFlow *gflow) : Base(gflow) {};

    //! \brief Do whatever corrections the body particles need.
    virtual void correct() = 0;
  };

}
#endif // __BODY_HPP__GFLOW__