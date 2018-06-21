#ifndef __SECTORIZATION_HPP__GFLOW__
#define __SECTORIZATION_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  class Sectorization : public Base {
  public:
    Sectorization(GFlow *);
    
    // GFlow is a friend class
    friend class GFlow;
    
  private:
    
  };

};

#endif // __SECTORIZATION_HPP__GFLOW__