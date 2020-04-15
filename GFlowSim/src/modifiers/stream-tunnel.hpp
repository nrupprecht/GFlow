#ifndef __STREAM_TUNNEL_HPP__GFLOW__
#define __STREAM_TUNNEL_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../utility/randomengines.hpp"

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
    //! \brief Function that is called recursively to fill the entry area with particles, no matter the dimensionality of the simulation.
    inline void recursive_fill(const int, Vec&, Vec&) const;

    //! \brief The amount of space in which we create and push particles.
    real entry_width = 2.f;
    //! \brief The amount of space in which we stablilize particle velocity at the end.
    real exit_width = 2.f;
    //! \brief The fraction of the entry width that should be used to add particles.
    real entry_fraction = 0.25f;

    //! \brief The velocity at which the particles should be driven.
    real driving_velocity = 1.f;
    //! \brief A vec that stores the velocity, { driving_velocity, 0, 0, ... }.
    Vec driving_velocity_vec;

    //! The min radius of an added particle.
    real min_r = 0.05f;
    //! The max radius of an added particle.
    real max_r = 0.05f;
    //! \brief The target system (volume) density.
    real phi_target = 0.9f;
    //! \brief The spacing needed between adjacent particles.
    real ave_spacing = 0.05f;

    //! \brief A random engine for creating particles.
    //!
    //! Because of C++, and how the constructor is set up, this needs to be after min_r and max_r in the class definition.
    ProportionalRandomEngine random_radius;

    // --- Internal parameters ---

    //! \brief The place at which particles should start being effected by the velocity field.
    real entry_threshold;
    //! \brief The place at which particles should stop being effected by the velocity field.
    real exit_threshold;

    //! \brief The amount of time since we last created particles.
    real last_creation_time = 0.f;

    //! \brief The x position where we should start placing particles.
    real next_x_coord;
  };

}
#endif // __STREAM_TUNNEL_HPP__GFLOW__