#ifndef __PARTICLE_FIXER_HPP__GFLOW__
#define __PARTICLE_FIXER_HPP__GFLOW__

#include "../utility/randomengines.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief An object that keeps track of how to initialize the velocity of a particle.
  *
  */
  struct ParticleFixer {
    //! \brief Constructor.
    ParticleFixer(int d, int id) : global_id(id), dimensions(d) {
      velocity = new RealType[d];
    };

    ParticleFixer(const ParticleFixer &fix) : global_id(fix.global_id), dimensions(fix.dimensions) {
      velocity = new RealType[dimensions];
      copyVec(fix.velocity, velocity, dimensions);
    }

    ParticleFixer(ParticleFixer &&fix) : global_id(fix.global_id), dimensions(fix.dimensions) {
      velocity = fix.velocity;
      fix.velocity = nullptr;
    }

    //! \brief Destructor.
    ~ParticleFixer() {
      if (velocity) delete [] velocity;
    }

    ParticleFixer& operator=(const ParticleFixer &fix) {
      global_id = fix.global_id;
      // Handle velocity vector
      if (dimensions!=fix.dimensions) {
        if (velocity) delete [] velocity;
        velocity = new RealType[dimensions];
      }
      copyVec(fix.velocity, velocity, dimensions);
      // Handle dimensions
      dimensions = fix.dimensions;
      // Return 
      return *this;
    }

    ParticleFixer& operator=(ParticleFixer &&fix) {
      global_id = fix.global_id;
      // Handle velocity vector
      velocity = fix.velocity;
      fix.velocity = nullptr;
      // Handle dimensions
      dimensions = fix.dimensions;
      // Return 
      return *this;
    }

    //! \brief The initial velocity.
    RealType *velocity = nullptr;

    //! \brief The global id of the particle that this object will fix.
    int global_id;

    //! \brief The dimensionality of the particle fixer.
    int dimensions;
  };

}
#endif // __PARTICLE_FIXER_HPP__GFLOW__