#ifndef __PARTICLE_FIXER_HPP__GFLOW__
#define __PARTICLE_FIXER_HPP__GFLOW__

#include "../utility/randomengines.hpp"
#include "../utility/vec.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief An object that keeps track of how to initialize the velocity of a particle.
  *
  */
  struct ParticleFixer {
    //! \brief Constructor.
    ParticleFixer(int d, int id) : global_id(id), velocity(d) {};

    ParticleFixer(const ParticleFixer &fix) : global_id(fix.global_id), velocity(fix.velocity) {};

    ParticleFixer(ParticleFixer &&fix) : global_id(fix.global_id), velocity(0) {
      velocity = std::move(fix.velocity);
      fix.velocity.data = nullptr;
    }

    ParticleFixer& operator=(const ParticleFixer &fix) {
      global_id = fix.global_id;
      // Handle velocity vector
      velocity = fix.velocity;
      // Return 
      return *this;
    }

    ParticleFixer& operator=(ParticleFixer &&fix) {
      global_id = fix.global_id;
      // Handle velocity vector
      velocity = std::move(fix.velocity);
      fix.velocity.data = nullptr;
      // Return 
      return *this;
    }

    //! \brief The initial velocity
    Vec velocity;

    //! \brief The global id of the particle that this object will fix.
    int global_id;
  };

}
#endif // __PARTICLE_FIXER_HPP__GFLOW__