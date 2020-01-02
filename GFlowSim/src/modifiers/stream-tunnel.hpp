#ifndef __STREAM_TUNNEL_HPP__GFLOW__
#define __STREAM_TUNNEL_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {
  
  class StreamTunnel : public Modifier {
  public:

    //! \brief Enforce wind tunnel conditions.
    virtual void post_forces() override;

  private:

  };

}
#endif // __STREAM_TUNNEL_HPP__GFLOW__