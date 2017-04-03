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

// Wall constants
const double default_wall_repulsion     = 100.;
const double default_wall_dissipation   = 50.;
const double default_wall_coeff         = sqrt(0.5);

// Drag for findPackedSolution
const double default_packed_drag        = 5.;

#endif 
