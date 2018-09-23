#ifndef __MODIFIER_HPP__
#define __MODIFIER_HPP__

#include "../gflow.hpp"
#include "simdata.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  class Modifier : public Base {
  public:
    // Constructor
    Modifier(GFlow *);

    bool getRemove();

    void setRemove(bool);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    bool remove;
  };

}
#endif // __MODIFIER_HPP__