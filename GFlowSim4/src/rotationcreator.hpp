#ifndef __ROTATION_CREATOR_HPP__GFLOW__
#define __ROTATION_CREATOR_HPP__GFLOW__

#include "creator.hpp"

namespace GFlowSimulation {

  class RotationCreator : public Creator {
  public:
    // Constructor -- pass in command line arguments
    RotationCreator(int, char**);

    // Constructor -- pass in an ArgParse
    RotationCreator(ArgParse *parser);

    virtual class GFlow* createSimulation() final;
  };

}
#endif // __ROTATION_CREATOR_HPP__GFLOW__