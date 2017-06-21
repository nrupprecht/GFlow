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
    // Constructors
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

  /*
   * @class Insertion
   *
   */
  class Insertion : public Characteristic {
  public:
    // Constructors
    Insertion(vec2 v) : velocity(v), omega(0), distance(0), speed(sqrt(sqr(v))), forward(true), useV(true), useOm(false) {};
    Insertion(vec2 v, RealType om) : velocity(v), omega(om), distance(0), speed(sqrt(sqr(v))), forward(true), useV(true), useOm(false) {};
    
    virtual void modify(SimData*, int, RealType);
  private:
    vec2 velocity;
    RealType omega;
    RealType distance, speed;
    bool forward;
    bool useV, useOm;
  };

}
#endif // __CHARACTERISTICS_HPP__
