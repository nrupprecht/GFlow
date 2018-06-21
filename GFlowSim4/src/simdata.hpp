#ifndef __SIM_DATA_HPP__GFLOW__
#define __SIM_DATA_HPP__GFLOW__

#include "gflow.hpp"
#include "utility.hpp"

namespace GFlowSimulation {

  /*
  *  @class SimData
  *
  *  Contains the data for the particles in the simulation
  *
  */
  class SimData : protected Base {
  public:
    // Constructor
    SimData(GFlow *);

    // Destructor
    ~SimData();

    // Clean all pointers
    void clean();

    // GFlow is a friend class
    friend class GFlow;

    // Particle data

    // Number of particles
    int number;

    // Number of "ghost" particles
    int numberG;

    // All particles have position, force - this data is used for integration, calculating verlet lists, etc.
    RealType **x, **v, **f;
    // All particles also have a characteristic length (sg - for sigma), (inverse) mass
    RealType *sg, *im;
    // Store the type of each atom
    int *type;
    // More data, for more complex objects
    RealType **dataF; // Floating point data
    RealType **dataI; // Integer data

    // Periodicity (wrapping) indicator
    bool wrap[DIMENSIONS];
  };


}
#endif // __SIM_DATA_HPP__GFLOW__