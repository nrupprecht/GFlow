/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 12, 2017
 *
 */

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
  const RealType default_wall_width           = 0.0;

  // Sectorization constants
  const RealType default_sectorization_skin_depth = 0.025;

  // Integrator constants
  const RealType default_epsilon              = 1e-4;
  const RealType default_max_timestep         = 5e-3;
  const RealType default_min_timestep         = 5e-6;
  const RealType default_update_delay         = 0.002;
  const RealType default_max_update_delay     = 0.01;
  const RealType default_delay_factor         = 0.9;
  const int default_period_iterations         = 150;

  // Scalar field constants
  const double default_diffusion          = 1.;
  const double default_lambda             = 1.;

}
#endif // __DEFAULT_CONSTANTS_HPP__
