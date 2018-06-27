#ifndef __INTEGRATOR_HPP__GFLOW__
#define __INTEGRATOR_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  class Integrator : public Base {
  public:
    // Constructor
    Integrator(GFlow *);
    
    // --- Accessors
    RealType getTimeStep();

    // --- Mutators
    void setDT(RealType);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // Time step
    RealType dt;

  };

}
#endif // __INTEGRATOR_HPP__GFLOW__