#ifndef DEFAULT_CONSTANTS_H
#define DEFAULT_CONSTANTS_H

#include "Utility.h"

// Sphere (Disk) constants
const double default_sphere_density     = 1.;
const double default_sphere_repulsion   = 100.;
const double default_sphere_dissipation = 7.5;  // If this is to big (compaired to epsilon, presumably), BAD things happen
const double default_sphere_coeff       = sqrt(0.3);

// Lenard Jones constants
const double LJ_cutoff_factor           = 2.5;

// Electric field constants
const double EF_cutoff_factor           = 0.1 ; // Distance per unit charge

// Wall constants
const double default_wall_repulsion     = 100.;
const double default_wall_dissipation   = 50.;
const double default_wall_coeff         = sqrt(0.5);

// Drag for findPackedSolution
const double default_packed_drag        = 5.;

// Bacteria constants
const double default_bacteria_reorient = 0.1;
const double default_bacteria_strength = 0.1;
const double default_bacteria_drag     = 10.;
const double default_bacteria_target_velocity = 0.25;
const double default_bacteria_target_velocity_sqr = sqr(default_bacteria_target_velocity);
const double default_bacteria_reproduction_const = 5.;
const double default_bacteria_death_const = 10.;

// Other
const double default_upper_window_factor = 6;

#endif 
