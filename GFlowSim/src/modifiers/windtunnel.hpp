#ifndef __WIND_TUNNEL_HPP__GFLOW__
#define __WIND_TUNNEL_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class WindTunnel : public Modifier {
  public:
    WindTunnel(GFlow*, RealType);

    virtual void post_forces() override;


  private:

    //! @brief The place at which particles should stop being effected by the velocity field.
    RealType leftBound;
    //! @brief The place at which particles should start being effected by the velocity field.
    RealType rightBound;

    //! @brief The velocity given to particles that are transfered.
    RealType velocity;
  };

}
#endif // __WIND_TUNNEL_HPP__GFLOW__