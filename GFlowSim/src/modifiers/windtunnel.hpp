#ifndef __WIND_TUNNEL_HPP__GFLOW__
#define __WIND_TUNNEL_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  //! \brief Simulates a ``wind tunnel'' by forcing particles near the left and right edges of the simulation
  //! to move at a specific velocity (or at least very near it). This drives particle flow through the simulation.
  class WindTunnel : public Modifier {
  public:
    //! \brief Constructor.
    WindTunnel(GFlow*, RealType);
    //! \brief Constructor.
    WindTunnel(GFlow*);

    //! \brief Enforce wind tunnel conditions.
    virtual void post_forces() override;

    //! \brief Create this object from parse node data.
    virtual void parse_construct(HeadNode*, const std::map<string, string>&) override;

  private:
    //! \brief The half-width of the area in which particles are corralled by the velocity field.
    RealType halfWidth = 4.f;

    //! \brief The place at which particles should stop being effected by the velocity field.
    RealType leftBound = 0.f;

    //! \brief The place at which particles should start being effected by the velocity field.
    RealType rightBound = 0.f;

    //! \brief The target acceleration for the particles in the velocity field.
    RealType acceleration = 2.f;

    //! \brief The velocity given to particles that are transfered.
    RealType velocity = 2.f;
  };

}
#endif // __WIND_TUNNEL_HPP__GFLOW__