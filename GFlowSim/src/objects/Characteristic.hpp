#ifndef __CHARACTERISTICS_HPP__
#define __CHARACTERISTICS_HPP__

// Includes
#include "../../include/vec2d.hpp"

namespace GFlow {

  // Forward declaration to SimData
  class SimData;
  

  /*
   * @class Characteristic
   *
   */
  class Characteristic {
  public:
    virtual void modify(SimData*, int, RealType) = 0;
  };

  /*
   * @class ConstantVelocity
   *
   */
  class ConstantVelocity : public Characteristic {
  public:
    ConstantVelocity(vec2 v) : velocity(v), omega(0), useV(true), useOm(false), active(true) {};
    ConstantVelocity(RealType om) : velocity(Zero), omega(om), useV(false), useOm(true), active(true) {};
    ConstantVelocity(vec2 v, RealType om) : velocity(v), omega(om), useV(true), useOm(true), active(true) {};

    virtual void modify(SimData*, int, RealType);
  private:
    vec2 velocity;
    RealType omega;
    bool useV, useOm;

    bool active;
  };

}
#endif // __CHARACTERISTICS_HPP__
