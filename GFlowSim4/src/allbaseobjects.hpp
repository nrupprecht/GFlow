#ifndef __ALL_BASE_OBJECTS_HPP__GFLOW__
#define __ALL_BASE_OBJECTS_HPP__GFLOW__

#include "simdata.hpp"
// Base class for integrators
#include "integrator.hpp"
// Base class for interactions/forces
#include "interaction.hpp"
// Specific domains
#include "sectorization.hpp"
#include "domain.hpp"

#include "domain2d_test.hpp" // A test

// Other ...
#include "verletlist.hpp"
#include "communicator.hpp"
#include "datamaster.hpp"
#include "forcemaster.hpp"
#include "modifier.hpp"

#endif