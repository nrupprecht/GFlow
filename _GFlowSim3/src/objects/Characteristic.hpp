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
    // Default constructor
    Characteristic() : active(true) {};
    // Characteristic modifies particle or wall, bool is false for walls
    virtual void modify(SimData*, int, RealType, bool=true) = 0;
    // Reset the characteristic if neccessary
    virtual void reset() {};

    // Get active flag
    bool getActive() { return active; }

    // Set active flag
    void setActive(bool a) { active = a; }

  protected:
    // Whether the characteristic is active or not
    bool active;
  };

  /*
   * @class ConstantVelocity
   *
   */
  class ConstantVelocity : public Characteristic {
  public:
    // Constructors
    ConstantVelocity(vec2 p, vec2 v) : pos(p), time(0), velocity(v), omega(0), useV(true), useOm(false), stop(false) {};
    ConstantVelocity(vec2 p, vec2 v, bool s) : pos(p), time(0), velocity(v), omega(0), useV(true), useOm(false), stop(s) {};
    ConstantVelocity(RealType om) : pos(Zero), time(0), velocity(Zero), omega(om), useV(false), useOm(true), stop(true) {};
    ConstantVelocity(vec2 p, vec2 v, RealType om) : pos(p), time(0), velocity(v), omega(om), useV(true), useOm(true), stop(false) {};
    ConstantVelocity(vec2 p, vec2 v, RealType om, bool s) : pos(p), time(0), velocity(v), omega(om), useV(true), useOm(true), stop(s) {};

    virtual void modify(SimData*, int, RealType, bool=true);
    virtual void reset();
  private:

    vec2 pos; // Initial position
    RealType time; // How much time has passed

    vec2 velocity;
    RealType omega;
    // Whether to control velocity, angular velocity
    bool useV, useOm;
    
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
    
    virtual void modify(SimData*, int, RealType, bool=true);
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

    virtual void modify(SimData*, int, RealType, bool=true);
  private:
    RealType radius, omega, theta;
    vec2 center;
    bool first;
  };

  class Fixed : public Characteristic {
  public:
    // Constructors
    Fixed(vec2 p) : pos(p) {};
    Fixed(int, SimData*, bool);
    
    virtual void modify(SimData*, int, RealType, bool=true);
  private:
    vec2 pos; // The position the particle should stay in
  };

  class ApplyForce : public Characteristic {
  public:
    // Constructors
    ApplyForce(vec2);
    ApplyForce(vec2, vec2);
    
    virtual void modify(SimData*, int, RealType, bool=true);
    // Reset the characteristic if neccessary
    virtual void reset();
  private:
      vec2 F0, F, dF;
  };

}
#endif // __CHARACTERISTICS_HPP__
