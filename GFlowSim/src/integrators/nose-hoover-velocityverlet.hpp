#ifndef __NOSE_HOOVER_VELOCITY_VERLET_HPP__GFLOW__
#define __NOSE_HOOVER_VELOCITY_VERLET_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  template<int dimensions> class NoseHooverVelocityVerlet : public Integrator {
  public:
    // Constructor
    NoseHooverVelocityVerlet(GFlow *);
    virtual void pre_forces() override;
    virtual void post_forces() override;

  private:
    inline void update_positions();
    inline void update_velocities();

    //! \brief The target value of T.
    RealType temperature;

    //! \brief The Nose-Hoover zeta value
    RealType noseHooverZeta = 0;

    //! \brief The Nose-Hoover relaxation factor.
    RealType noseHooverQ;

    //! \brief The amount of time to wait before assigning a target temperature, if one has not already been set.
    RealType wait_time;
  };

  // Include the implementation file.
  #include "nose-hoover-velocity-verlet.tpp"

  inline Integrator* choose_nose_hoover_velocity_verlet(GFlow *gflow, int sim_dimensions) {
    switch (sim_dimensions) {
      case 1:
        return new NoseHooverVelocityVerlet<1>(gflow);
      case 2:
        return new NoseHooverVelocityVerlet<2>(gflow);
      case 3: 
        return new NoseHooverVelocityVerlet<3>(gflow);
      case 4: 
        return new NoseHooverVelocityVerlet<4>(gflow);
      default:
        throw BadDimension();
    }
  }

}
#endif // __NOSE_HOOVER_VELOCITY_VERLET_HPP__GFLOW__