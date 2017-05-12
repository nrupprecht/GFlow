#ifndef __DEFAULT_CONSTANTS_HPP__
#define __DEFAULT_CONSTANTS_HPP__

#include "Utility.hpp"

namespace GFlow {

  // Particle constants
  const RealType default_particle_density     = 1.;
  const RealType default_particle_repulsion   = 150.;
  const RealType default_particle_dissipation = 15.;
  const RealType default_particle_coeff       = sqrt(0.3); 
  
  // Wall constants
  const RealType default_wall_repulsion       = 150.;
  const RealType default_wall_dissipation     = 50.;
  const RealType default_wall_coeff           = sqrt(0.5);

}
#endif // __DEFAULT_CONSTANTS_HPP__
