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
  ParticleFixer(int d, int id)
      : velocity(d), global_id(id) {};

  //! \brief The initial velocity
  Vec velocity;

  //! \brief The global id of the particle that this object will fix.
  int global_id;
};

}
#endif // __PARTICLE_FIXER_HPP__GFLOW__