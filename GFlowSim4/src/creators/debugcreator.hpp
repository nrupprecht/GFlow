#ifndef __DEBUG_CREATOR_HPP__GFLOW__
#define __DEBUG_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"

namespace GFlowSimulation {

  class DebugCreator : public Creator {
  public:
    // Constructor
    DebugCreator(int, char**);

    // Constructor
    DebugCreator(ArgParse*);

    // Create simulation
    virtual GFlow* createSimulation();
  };

}
#endif // __DEBUG_CREATOR_HPP__GFLOW__