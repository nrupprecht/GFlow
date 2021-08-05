#ifndef __CREATOR_HPP__GFLOW__
#define __CREATOR_HPP__GFLOW__

#include "../utility/utility.hpp"
#include "../utility/ArgParse.hpp"
#include "../utility/vectormath.hpp"
#include "../gflow.hpp"
#include "../creators/particlefixer.hpp"

// Accumulation files
#include "../allbaseobjects.hpp"
#include "../alldataobjects.hpp"
#include "../allintegrators.hpp"
#include "../allmodifiers.hpp"
#include "../allinteractions.hpp"

namespace GFlowSimulation {

/*
*  \brief Base class for simulation creators
*
*  Creates GFlow objects that can then be run. All creator objects inherit
*  from this base class.
*/
class Creator {
 public:
  //! \brief Constructor -- pass in command line arguments.
  Creator(int, char **);

  //! \brief Constructor -- pass in a pointer to an ArgParse object.
  Creator(ArgParse *);

  //! \brief Destructor.
  virtual ~Creator();

  //! \brief Set the random seed.
  void setSeed(uint);

  //! \brief Get the random seed.
  unsigned getSeed();

  //! \brief Seed whatever random generators there are.
  virtual void seedGenerator(uint);

  //! \brief Create a GFlow Object
  virtual class GFlow *createSimulation() = 0;

  //! \brief Set the boundary condition flags.
  void setBCFlag(BCFlag b);

  //! \brief Set the dimensionality of the simulation being created.
  void setDimensions(int);

  //! \brief Hard sphere version of the relaxation function.
  //!
  //! This function gives the simdata an overdamped integrator and hard sphere interactions, runs it for some amount
  //! of time, then discards the integrator and resets the timers.
  static void hs_relax(class GFlow *, RealType= 0.25, bool= true, HeadNode * = nullptr);

  //! \brief Native version of the relaxation function.
  //!
  //! Relaxation function that uses the native forces and an overdamped integrator. It runs for some amount of time,
  //! the discards the integrator and resets the timers.
  static void relax(class GFlow *, RealType= 0.25, HeadNode * = nullptr);

  //! \brief Only needed for multiprocessor MPI runs. Corrects the global ids of particles on different processors so that they are all unique.
  static void correct_global_ids(class GFlow *);

  //! \brief Use the particle fixers to assign particle velocities.
  void fix_particle_velocities(const shared_ptr<SimData> &);

  //! \brief Clear all the particle fixers
  void clear_particle_fixers();

 protected:

  // Command line arguments
  //! \brief The number of command line arguments. Has to be set from the outside.
  int argc;
  //! \brief The command line arguments. Has to be set from the outside.
  char **argv;

  // --- Data that all creators will (likely) use

  //! \brief Boundary condition flag.
  //!
  //! Use for when we want all the boundary conditions to be the same.
  BCFlag bcFlag = BCFlag::WRAP;

  //! \brief Bounds of the simulation we are creating.
  Bounds simBounds;

  //! \brief Argument parser.
  ArgParse *parserPtr;

  //! \brief Parser flag.
  //!
  //! True if we allocated the parser ourselves. False if a parser was passed in to us.
  bool ourParser;

  //! \brief The random seed.
  unsigned seed;

  //! \brief The dimensionality of the simulation.
  int sim_dimensions;

  //! \brief Vector of particle fixers. These are used to assign velocities to particles after the relax step.
  mutable vector<ParticleFixer> particle_fixers;
};

}
#endif // __CREATOR_HPP__GFLOW__
