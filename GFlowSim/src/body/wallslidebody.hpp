#ifndef __WALL_SLIDE_BODY_HPP__GFLOW__
#define __WALL_SLIDE_BODY_HPP__GFLOW__

#include "../base/body.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class WallSlideBody : public Body, public Group {
  public:
    //! \brief Default constructor. Can set the slide dimension of the body.
    WallSlideBody(GFlow*, int=0);

    //! \brief Constructor that takes a group of particles.
    WallSlideBody(GFlow*, Group, int=0);

    //! \brief Set all the velocities of the group particles to zero.
    virtual void pre_integrate() override;

    //! \brief Do whatever corrections the body particles need.
    virtual void correct() override;

  private:
    //! \brief What dimension the body is allowed to slide in.
    int slide_dimension = 0;
  };

}
#endif // __WALL_SLIDE_BODY_HPP__GFLOW__