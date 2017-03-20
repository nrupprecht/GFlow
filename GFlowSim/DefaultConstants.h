#ifndef DEFAULT_CONSTANTS_H
#define DEFAULT_CONSTANTS_H

// Sphere (Disk) constants
const double default_sphere_density     = 1.;
const double default_sphere_repulsion   = 50.;
const double default_sphere_dissipation = 10.;
const double default_sphere_coeff       = sqrt(0.5);
const double default_sphere_drag        = 1.;

// Wall constants
const double default_wall_repulsion     = 100.;
const double default_wall_dissipation   = 50.;
const double default_wall_coeff         = sqrt(0.5);

// Drag for findPackedSolution
const double default_packed_drag        = 5.;

typedef float floatType; // What type of floating point should we use for our data
typedef vect<floatType> vec2;

#endif 
