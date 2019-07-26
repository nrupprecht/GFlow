#ifndef __ALL_FORCES_HPP__GFLOW__
#define __ALL_FORCES_HPP__GFLOW__

// Parent classes
#include "interactions/hard_sphere.hpp"
#include "interactions/hard_sphere_ds.hpp"
#include "interactions/hard_sphere_cf.hpp"
#include "interactions/lennard_jones.hpp"
#include "interactions/buckingham.hpp"

// Specific classes
#include "interactions/hard_sphere__2d.hpp"
#include "interactions/hard_sphere__3d.hpp"

#include "interactions/hard_sphere_ds__2d.hpp"
#include "interactions/hard_sphere_ds__3d.hpp"

#include "interactions/hard_sphere_cf__2d.hpp"

#include "interactions/hard_sphere__reflecting__2d.hpp"

#include "interactions/lennard_jones__2d.hpp"
#include "interactions/lennard_jones__3d.hpp"

#include "interactions/coulomb.hpp"

#include "interactions/buckingham__2d.hpp"
#include "interactions/buckingham__3d.hpp"

#include "interactions/detector.hpp"
#include "interactions/detector__2d.hpp"
#include "interactions/detector__3d.hpp"

// Specialized
#include "interactions/demon_wall.hpp"

#endif // __ALL_FORCES_HPP__GFLOW__