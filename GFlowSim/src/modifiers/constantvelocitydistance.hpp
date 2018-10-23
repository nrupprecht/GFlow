#ifndef __CONSTANT_VELOCITY_DISTANCE__GFLOW__
#define __CONSTANT_VELOCITY_DISTANCE__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  /*
  *  @class ConstantVelocity
  *
  *  A modifier that forces an object to move at constant velocity.
  */
  class ConstantVelocityDistance : public Modifier {
  public:
    //! @brief Constructor
    ConstantVelocityDistance(GFlow*, int, RealType*, RealType);

    //! @brief Set force to zero and velocity to the specified velocity.
    virtual void post_forces() override;

  private:
    //! @brief The global id of the particle this modifier modifies.
    int global_id;

    //! @brief The constant velocity.
    RealType velocity[DIMENSIONS];

    //! @brief The displacement covered so far.
    RealType displacement[DIMENSIONS];

    //! @brief The distance to cover before stopping.
    RealType distance;

    //! @brief Whether we still have to move.
    bool moving = true;
  };

}
#endif // __CONSTANT_VELOCITY_DISTANCE__GFLOW__