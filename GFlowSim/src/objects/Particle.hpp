#ifndef __PARTICLE_HPP__
#define __PARTICLE_HPP__

// Includes
#include "../../include/vec2d.hpp"
#include "../../include/CSVUtility.hpp"
#include "../../include/DefaultConstants.hpp"

namespace GFlow {
  
  /*
   * @class Particle
   *
   */
  struct Particle {
    // Default constructor
    Particle();

    // Initialize position with vector and sigma
    Particle(vec2, RealType);

    // Initialize position with reals and sigma
    Particle(RealType, RealType, RealType);

    vec2 position, velocity, force;
    RealType theta, omega, torque;

    RealType sigma, invMass, invII;
    RealType repulsion, dissipation, coeff;

    int interaction;
  };

  inline string toCSV(const Particle& p) {
    stringstream stream;
    string str;
    stream << p.position.x << "," << p.position.y << "," << p.sigma << "," << p.theta << "," << p.interaction << ",0";
    stream >> str;
    return str;
  }
  
}

#endif // __PARTICLE_HPP__
