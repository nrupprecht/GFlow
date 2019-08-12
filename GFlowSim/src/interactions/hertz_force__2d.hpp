#ifndef __HERTZ_FORCE__2D_HPP__GFLOW__
#define __HERTZ_FORCE__2D_HPP__GFLOW__

#include "hertz_force.hpp"
#include "interaction__2d.hpp"

namespace GFlowSimulation {

  class HertzForce2d : public HertzForce, public Interaction2d {
  public:

  private:
    void kernel(int, int, RealType, RealType, RealType, RealType*, RealType**) const override;

  };
}
#endif // __HERTZ_FORCE__2D_HPP__GFLOW__