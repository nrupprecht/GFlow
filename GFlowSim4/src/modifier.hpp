#ifndef __MODIFIER_HPP__
#define __MODIFIER_HPP__

#include "gflow.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  class Modifier : public Base {
  public:
    // Constructor
    Modifier(GFlow *);

    // GFlow is a friend class
    friend class GFlow;
  };

}
#endif // __MODIFIER_HPP__