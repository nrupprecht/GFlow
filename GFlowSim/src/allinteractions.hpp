#ifndef __ALL_FORCES_HPP__GFLOW__
#define __ALL_FORCES_HPP__GFLOW__

// Parent classes
#include "interactions/hard_sphere.hpp"
#include "interactions/hard_sphere_ds.hpp"
#include "interactions/lennard_jones.hpp"
#include "interactions/buckingham.hpp"

// Specific classes
#include "interactions/hard_sphere__verlet_pairs__2d.hpp"
#include "interactions/hard_sphere__verlet_pairs__3d.hpp"

#include "interactions/hard_sphere_ds__verlet_pairs__2d.hpp"
#include "interactions/hard_sphere_ds__verlet_pairs__3d.hpp"

#include "interactions/lennard_jones__verlet_pairs__2d.hpp"
#include "interactions/lennard_jones__verlet_pairs__3d.hpp"

#include "interactions/buckingham__verlet_pairs__2d.hpp"
#include "interactions/buckingham__verlet_pairs__3d.hpp"

#endif // __ALL_FORCES_HPP__GFLOW__