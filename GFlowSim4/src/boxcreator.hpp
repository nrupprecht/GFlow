#ifndef __BOX_CREATOR_HPP__GFLOW__
#define __BOX_CREATOR_HPP__GFLOW__

#include "creator.hpp"

namespace GFlowSimulation {

  class BoxCreator : public Creator {
  public:
    // Constructor
    BoxCreator(int, char**);

    // Constructor
    BoxCreator(ArgParse*);

    // Create simulation
    virtual GFlow* createSimulation();
  };

}
#endif