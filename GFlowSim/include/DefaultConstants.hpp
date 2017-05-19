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

  // Sectorization constants
  const RealType default_sectorization_skin_depth = 0.025;

  // Integrator constants
  const RealType default_epsilon              = 1e-4;
  const RealType default_update_delay         = 0.002;
  const RealType default_max_update_delay     = 0.01;
  const RealType default_delay_factor         = 1.3;
  const int default_period_iterations         = 150;
  const RealType default_max_timestep         = 5e-3;

}
#endif // __DEFAULT_CONSTANTS_HPP__
