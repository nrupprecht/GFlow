#ifndef __CONTAINER_LAYOUT_HPP__GFLOW__
#define __CONTAINER_LAYOUT_HPP__GFLOW__

// Include the relevant files.
#include "particle-data-soa.hpp"
#include "particle-data-aos.hpp"

namespace GFlowSimulation {
  // Choose which type of particle container to use.
  template<int dims> using ParticleContainer = ParticleContainer_SOA<dims>;
}

#endif // __CONTAINER_LAYOUT_HPP__GFLOW__