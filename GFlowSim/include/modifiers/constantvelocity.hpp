#ifndef __CONSTANT_VELOCITY_HPP__GFLOW__
#define __CONSTANT_VELOCITY_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  /*
  *  @class ConstantVelocity
  *
  *  A modifier that forces an object to move at constant velocity.
  */
  class ConstantVelocity : public Modifier {
  public:
    //! @brief Constructor.
    ConstantVelocity(GFlow*, int, RealType*);

    //! @brief Destructor.
    ~ConstantVelocity();

    //! @brief Set force to zero and velocity to the specified velocity.
    virtual void post_forces() override;

  private:
    //! @brief The global id of the particle this modifier modifies.
    int global_id;

    //! @brief The constant velocity.
    RealType *velocity;
  };

}
#endif // __CONSTANT_VELOCITY_HPP__GFLOW__