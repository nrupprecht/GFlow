#ifndef __CREATOR_HPP__GFLOW__
#define __CREATOR_HPP__GFLOW__

#include "simdata.hpp"
#include "integrator.hpp"
#include "sectorization.hpp"
#include "communicator.hpp"
#include "datamaster.hpp"
#include "modifier.hpp"

// Accumulation files
#include "allforces.hpp"
#include "alldataobjects.hpp"
#include "allintegrators.hpp"

namespace GFlowSimulation {

  /*
  *  @class Creator
  *  Creates GFlow objects that can then be run. All creator objects inherit
  *  from this base class.
  *
  */
  class Creator {
  public:
    // Constructor -- pass in command line arguments
    Creator(int, char**);

    // Create a GFlow Object
    virtual GFlow* createSimulation() = 0;

  private:
    // Command line arguments
    int argc;
    char **argv;

    // --- Data that all creators will (likely) use

    // Bounds of the simulation to create
    Bounds simBounds;
  };

}
#endif