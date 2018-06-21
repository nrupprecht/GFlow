#include "boxcreator.hpp"

namespace GFlowSimulation {

  BoxCreator::BoxCreator(int argc, char **argv) : Creator(argc, argv) {};

  GFlow* BoxCreator::createSimulation() {
    // Create a new gflow object
    GFlow *gflow = new GFlow;

    // Set the bounds of the gflow object

    return gflow;
  }

}