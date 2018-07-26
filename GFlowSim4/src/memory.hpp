#ifndef __MEMORY_HPP__GFLOW__
#define __MEMORY_HPP__GFLOW__

namespace GFlowSimulation {

  class Memory {
  public:
    // Aligned memory allocation
    static template<typename T> T* aligned_alloc(uint, uint) {

    }
  };

}

#endif // __MEMORY_HPP__GFLOW__