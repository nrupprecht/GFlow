#ifndef __CREATOR_HPP__GFLOW__
#define __CREATOR_HPP__GFLOW__

#include "utility.hpp"
#include "ArgParse.hpp"
#include "vectormath.hpp"
#include "gflow.hpp"

// Accumulation files
#include "allbaseobjects.hpp"
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

    // Constructor -- pass in an ArgParse
    Creator(ArgParse *parser);

    // Destructor
    virtual ~Creator();

    // Set the random seed
    void setSeed(uint);

    // Get the random seed
    unsigned getSeed();

    // Seed whatever random generators there are
    virtual void seedGenerator(uint);

    // Create a GFlow Object
    virtual class GFlow* createSimulation() = 0;

    void setBCFlag(BCFlag b) { bcFlag = b; }

  protected:
    // --- Helper functions

    // Relax function - gives the simdata an overdamped integrator and hard sphere interactions, runs it for some amount
    // of time, then discards the integrator and resets the timers
    void hs_relax(class GFlow*, RealType=0.25);

    // Relax function - uses the native forces and an overdamped integrator
    void relax(class GFlow*, RealType=0.25);

    // Command line arguments
    int argc;
    char **argv;

    // --- Data that all creators will (likely) use

    // Boundary condition
    BCFlag bcFlag;

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