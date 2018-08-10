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
#include "allmodifiers.hpp"

namespace GFlowSimulation {

  /*
  *  @brief Base class for simulation creators
  *
  *  Creates GFlow objects that can then be run. All creator objects inherit
  *  from this base class.
  */
  class Creator {
  public:
    //! Constructor -- pass in command line arguments.
    Creator(int, char**);

    //! Constructor -- pass in a pointer to an ArgParse object.
    Creator(ArgParse*);

    //! Destructor.
    virtual ~Creator();

    //! Set the random seed.
    void setSeed(uint);

    //! Get the random seed.
    unsigned getSeed();

    //! Seed whatever random generators there are.
    virtual void seedGenerator(uint);

    //! Create a GFlow Object
    virtual class GFlow* createSimulation() = 0;

    //! Set the boundary condition flags.
    void setBCFlag(BCFlag b) { bcFlag = b; }

  protected:
    // --- Helper functions

    //! @brief Hard sphere version of the relaxation function.
    //!
    //! This function gives the simdata an overdamped integrator and hard sphere interactions, runs it for some amount
    //! of time, then discards the integrator and resets the timers.
    void hs_relax(class GFlow*, RealType=0.25);

    //! @brief Native version of the relaxation function.
    //!
    //! Relaxation function that uses the native forces and an overdamped integrator. It runs for some amount of time, 
    //! the discards the integrator and resets the timers.
    void relax(class GFlow*, RealType=0.25);

    // Command line arguments
    //! The number of command line arguments. Has to be set from the outside.
    int argc;
    //! The command line arguments. Has to be set from the outside.
    char **argv;

    // --- Data that all creators will (likely) use

    //! @brief Boundary condition flag.
    //!
    //! Use for when we want all the boundary conditions to be the same.
    BCFlag bcFlag;

    //! Bounds of the simulation we are creating.
    Bounds simBounds;

    //! Argument parser.
    ArgParse *parserPtr;

    //! @brief Parser flag.
    //!
    //! True if we allocated the parser ourselves. False if a parser was passed in to us.
    bool ourParser; 

    //! The random seed.
    unsigned seed;
  };

}
#endif // __CREATOR_HPP__GFLOW__