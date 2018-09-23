#include "rotationcreator.hpp"

namespace GFlowSimulation {

  RotationCreator::RotationCreator(int argc, char **argv) : Creator(argc, argv) {};

  RotationCreator::RotationCreator(ArgParse *parser) : Creator(parser) {};

  GFlow* RotationCreator::createSimulation() {
    return nullptr;
  }

}