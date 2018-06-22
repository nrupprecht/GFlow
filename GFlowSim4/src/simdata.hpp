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
  class SimData : public Base {
  public:
    // Constructor
    SimData(GFlow *);

    // Destructor
    ~SimData();

    // Clean all pointers
    void clean();

    // --- Functions for managing a SimData's data

    // Destroy old pointers, create arrays of a specified size
    void reserve(int);

    // Set all positions to zero
    void clearX();

    // Set all velocities to zero
    void clearV();

    // Set all forces to zero
    void clearF();

    // --- Particle data

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
    int      **dataI; // Integer data

    // The actual amount of data in the arrays
    int size;

    // GFlow is a friend class -- this allows it to access protected Base members
    friend class GFlow;
  };


}
#endif // __SIM_DATA_HPP__GFLOW__