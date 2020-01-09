#ifndef __STREAM_TUNNEL_HPP__GFLOW__
#define __STREAM_TUNNEL_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {
  
  class StreamTunnel : public Modifier {
  public:
    //! \brief Default constructor.
    StreamTunnel(GFlow*);
    //! \brief Constructor, sets velocity, min and max radii.
    StreamTunnel(GFlow*, const real, const real, const real);

    //! \brief Set the last creation time.
    virtual void pre_integrate() override;

    //! \brief Potentially add new particles.
    virtual void pre_forces() override;

    //! \brief Enforce wind tunnel conditions.
    virtual void post_forces() override;

    //! \brief Create this object from parse node data.
    virtual void parse_construct(HeadNode*, const std::map<string, string>&) override;

  private:
    //! \brief The amount of space in which we create and push particles.
    real entry_width = 4.f;
    //! \brief The amount of space in which we stablilize particle velocity at the end.
    real exit_width = 4.f;
    //! \brief The fraction of the entry width that should be used to add particles.
    real entry_fraction = 0.25f;

    //! \brief The place at which particles should start being effected by the velocity field.
    real entry_threshold;
    //! \brief The place at which particles should stop being effected by the velocity field.
    real exit_threshold;

    //! \brief The amount of time since we last created particles.
    real last_creation_time = 0.f;

    //! \brief The velocity at which the particles should be driven.
    real driving_velocity = 1.f;

    //! The min radius of an added particle.
    real min_r;
    //! The max radius of an added particle.
    real max_r;
    //! \brief The target system (volume) density.
    real phi_target = 0.9;
  };

}
#endif // __STREAM_TUNNEL_HPP__GFLOW__