#ifndef __MEMORY_OPTIMIZER_HPP__GFLOW__
#define __MEMORY_OPTIMIZER_HPP__GFLOW__

#include "utility.hpp"
#include "memory.hpp"

namespace GFlowSimulation {

  class MemoryOptimizer {
  public:
    static void GridParticles(class SimData&, const Bounds&);
  };

}
#endif // __MEMORY_OPTIMIZER_HPP__GFLOW__