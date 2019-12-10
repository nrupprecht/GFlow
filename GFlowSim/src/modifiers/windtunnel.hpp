#ifndef __WIND_TUNNEL_HPP__GFLOW__
#define __WIND_TUNNEL_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  //! \brief Simulates a ``wind tunnel'' by forcing particles near the left and right edges of the simulation
  //! to move at a specific velocity (or at least very near it). This drives particle flow through the simulation.
  class WindTunnel : public Modifier {
  public:
    //! @brief Constructor.
    WindTunnel(GFlow*, RealType);

    //! @brief Enforce wind tunnel conditions.
    virtual void post_forces() override;

  private:
    //! \brief The half-width of the area in which particles are corralled by the velocity field.
    RealType halfWidth;

    //! \brief The place at which particles should stop being effected by the velocity field.
    RealType leftBound;

    //! \brief The place at which particles should start being effected by the velocity field.
    RealType rightBound;

    //! \brief The target acceleration for the particles in the velocity field.
    RealType acceleration;

    //! \brief The velocity given to particles that are transfered.
    RealType velocity;
  };

}
#endif // __WIND_TUNNEL_HPP__GFLOW__