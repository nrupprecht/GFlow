#ifndef __WALL_SLIDE_BODY_HPP__GFLOW__
#define __WALL_SLIDE_BODY_HPP__GFLOW__

#include "base/body.hpp"
#include "other/group.hpp"

namespace GFlowSimulation {

class WallSlideBody : public Body, public Group {
 public:
  //! \brief Default constructor. Can set the slide dimension of the body.
  explicit WallSlideBody(GFlow* gflow, int= 0);

  //! \brief Constructor that takes a group of particles.
  WallSlideBody(GFlow* gflow, const Group& group, int d = 0);

  //! \brief Set all the velocities of the group particles to zero.
  virtual void pre_integrate() override;

  //! \brief Do whatever corrections the body particles need.
  virtual void correct() override;

  //! \brief Get the slide_dimension component of the wall position.
  RealType getPosition();

  //! \brief Get the slide_dimension component of the wall velocity.
  RealType getVelocity();

  //! \brief Get the length of the wall.
  RealType getLength() const;

  //! \brief Get the net force in the slide_dimension direction.
  RealType getFnet() const;

 private:
  //! \brief What dimension the body is allowed to slide in.
  int slide_dimension = 0;

  //! \brief The length of the wall.
  RealType length;

  //! \brief The net force in the slide_dimension direction.
  RealType Fnet = 0;
};

}
#endif // __WALL_SLIDE_BODY_HPP__GFLOW__