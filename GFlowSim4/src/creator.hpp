#ifndef __CREATOR_HPP__GFLOW__
#define __CREATOR_HPP__GFLOW__

#include "simdata.hpp"
#include "sectorization.hpp"
#include "communicator.hpp"
#include "datamaster.hpp"
#include "forcemaster.hpp"
#include "modifier.hpp"

// For relax step
#include "overdampedintegrator.hpp"

// Accumulation files
#include "allforces.hpp"
#include "alldataobjects.hpp"
#include "allintegrators.hpp"

// Argument parsing
#include "ArgParse.hpp"

// Useful
#include "vectormath.hpp"

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

    // Constructor -- pass in an ArgParse
    Creator(ArgParse *parser);

    // Destructor
    ~Creator();

    // Set the random seed
    void setSeed(uint);

    // Get the random seed
    unsigned getSeed();

    // Seed whatever random generators there are
    virtual void seedGenerator(uint);

    // Create a GFlow Object
    virtual GFlow* createSimulation() = 0;

  protected:
    // Command line arguments
    int argc;
    char **argv;

    // --- Data that all creators will (likely) use

    // Bounds of the simulation to create
    Bounds simBounds;

    // Parser
    ArgParse *parserPtr;
    bool ourParser; // True if we allocated the parser

    // Random seed
    unsigned seed;
  };

}
#endif