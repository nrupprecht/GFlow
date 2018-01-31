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
    // Characteristic modifies particle
    virtual void modify(SimData*, int, RealType) = 0;
    // Reset the characteristic if neccessary
    virtual void reset() {};
  };

  /*
   * @class ConstantVelocity
   *
   */
  class ConstantVelocity : public Characteristic {
  public:
    // Constructors
    ConstantVelocity(vec2 p, vec2 v) : pos(p), velocity(v), omega(0), useV(true), useOm(false), active(true), stop(false) {};
    ConstantVelocity(vec2 p, vec2 v, bool s) : pos(p), velocity(v), omega(0), useV(true), useOm(false), active(true), stop(s) {};
    ConstantVelocity(RealType om) : velocity(Zero), omega(om), useV(false), useOm(true), active(true), stop(true) {};
    ConstantVelocity(vec2 p, vec2 v, RealType om) : pos(p), velocity(v), omega(om), useV(true), useOm(true), active(true), stop(false) {};
    ConstantVelocity(vec2 p, vec2 v, RealType om, bool s) : pos(p), velocity(v), omega(om), useV(true), useOm(true), active(true), stop(s) {};

    virtual void modify(SimData*, int, RealType);
    virtual void reset();
  private:

    vec2 pos; // Initial position
    RealType time; // How much time has passed

    vec2 velocity;
    RealType omega;
    // Whether to control velocity, angular velocity
    bool useV, useOm;

    bool active;
    bool stop;
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

  /*
   * @class Circulate
   *
   */
  class Circulate : public Characteristic {
  public:
    // Constructors
    Circulate(RealType r, RealType om=2*PI, RealType th=0) : radius(r), omega(om), theta(r), center(Zero), first(true) {};

    virtual void modify(SimData*, int, RealType);
  private:
    RealType radius, omega, theta;
    vec2 center;
    bool first;
  };

  class Fixed : public Characteristic {
  public:
    // Constructors
    Fixed(vec2 p) : pos(p) {};
    
    virtual void modify(SimData*, int, RealType);
  private:
    vec2 pos; // The position the particle should stay in
  };

  class ApplyForce : public Characteristic {
  public:
    // Constructors
    ApplyForce(vec2 df) : F(Zero), dF(df) {};
    ApplyForce(vec2 f0, vec2 df) : F(f0), dF(df) {};
    
    virtual void modify(SimData*, int, RealType);
  private:
    vec2 F, dF;
  };

}
#endif // __CHARACTERISTICS_HPP__
