#ifndef __PARTICLE_HPP__
#define __PARTICLE_HPP__

namespace GFlow {
  
  /*
   * @class Particle
   *
   */
  struct Particle {
    Particle() {}; // STUBS

    vec2 position, velocity, force;
    RealType theta, omega, torque;

    RealType sigma, invMass, invII;
    RealType repulsion, dissipation, coeff;
  };
  
}

#endif __PARTICLE_HPP__
