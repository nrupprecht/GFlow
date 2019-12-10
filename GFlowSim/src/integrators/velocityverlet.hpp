#ifndef __VELOCITY_VERLET_HPP__GFLOW__
#define __VELOCITY_VERLET_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  template<int dimensions> class VelocityVerlet : public Integrator {
  public:
    // Constructor
    VelocityVerlet(GFlow *);
    virtual void pre_forces() override;
    virtual void post_forces() override;

  private:
    inline void update_positions();
    inline void update_velocities();
  };

  // Include the implementation file.
  #include "velocity-verlet.tpp"

  inline Integrator* choose_velocity_verlet(GFlow *gflow, int sim_dimensions) {
    switch (sim_dimensions) {
      case 1:
        return new VelocityVerlet<1>(gflow);
      case 2:
        return new VelocityVerlet<2>(gflow);
      case 3: 
        return new VelocityVerlet<3>(gflow);
      case 4: 
        return new VelocityVerlet<4>(gflow);
      default:
        throw BadDimension();
    }
  }

}
#endif // __VELOCITY_VERLET_HPP__GFLOW__