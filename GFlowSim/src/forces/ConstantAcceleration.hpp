#ifndef __CONSTANT_ACCELERATION_HPP__
#define __CONSTANT_ACCELERATION_HPP__

#include "ExternalForce.hpp"

namespace GFlow {
  
  /*
   * @class ConstantAcceleration
   * A constant acceleration, like gravity
   *
   */
  class ConstantAcceleration : public ExternalForce {
  public:
    // Default constructor - sets the viscosity to -9.81 m/s/s y
    ConstantAcceleration();

    // Acceleration setting constructor
    ConstantAcceleration(vec2);

  protected:
    // Inherited private virtual functions
    virtual void _applyForce(SimData*) const;
    virtual string _summary() const;
    // The acceleration vector
    vec2 acceleration;
  };

}

#endif // __CONSTANT_ACCELERATION_HPP__
