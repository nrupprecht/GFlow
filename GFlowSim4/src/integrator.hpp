#ifndef __INTEGRATOR_HPP__GFLOW__
#define __INTEGRATOR_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  class Integrator : protected Base {
  public:
    // Constructor
    Integrator(GFlow *);

    virtual void pre_step()    {};
    virtual void pre_forces()  {};
    virtual void post_forces() {};
    virtual void post_step()   {};

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // Time step
    RealType dt;

  };

}
#endif // __INTEGRATOR_HPP__GFLOW__