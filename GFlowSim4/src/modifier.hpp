#ifndef __MODIFIER_HPP__
#define __MODIFIER_HPP__

#include "gflow.hpp"

namespace GFlowSimulation {

  class Modifier : protected Base {
  public:
    // Constructor
    Modifier(GFlow *);

    virtual void pre_step()    {};
    virtual void pre_forces()  {};
    virtual void post_forces() {};
    virtual void post_step()   {};

    // GFlow is a friend class
    friend class GFlow;
  };

}
#endif // __MODIFIER_HPP__