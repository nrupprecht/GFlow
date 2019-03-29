#ifndef __RIGID_BODY_HPP__GFLOW__
#define __RIGID_BODY_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class RigidBody_2d : public Modifier, public Group {
  public:
    //! \brief Default constructor.
    RigidBody_2d(GFlow*);

    //! \brief Compute forces, torques, and update the particles' forces.
    virtual void post_forces() override;

    //! \brief Add the particle, and update the mass.
    virtual void add(int) override;

  protected:

    //! \brief The total mass of the body.
    RealType mass = 0;
  };

}
#endif // __RIGID_BODY_HPP__GFLOW__