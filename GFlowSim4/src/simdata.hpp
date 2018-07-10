#ifndef __SIM_DATA_HPP__GFLOW__
#define __SIM_DATA_HPP__GFLOW__

#include "gflow.hpp"
#include "utility.hpp"

namespace GFlowSimulation {

  /*
  *  @struct ParticleLayout
  *
  *  A way of encoding how the data should be laid out.
  *
  */
  struct ParticleLayout {
    ParticleLayout(string n, int t, int ndf, int ndi) : name(n), type(t), n_dataF(ndf), n_dataI(ndi) {};

    // --- Data
    string name;          // A name for the particle
    int type;             // What type this atom is denoted by
    int n_dataF, n_dataI; // The amount of dataF and dataI space allocated
  };

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

    // Destroy old pointers, create arrays of a specified size for x, v, f, sg, type
    void reserve(int);

  // Destroy old pointers, create arrays of a specified size for x, v, f, sg, type, dataF, dataI
    void reserveAll(int);

    // Set all positions to zero
    void clearX();

    // Set all velocities to zero
    void clearV();

    // Set all forces to zero
    void clearF();

    // --- Accessors

    RealType& X(int n, int d) { return x[n][d]; }
    RealType* X(int n)        { return x[n]; }
    RealType& V(int n, int d) { return v[n][d]; }
    RealType* V(int n)        { return v[n]; }
    RealType& F(int n, int d) { return f[n][d]; }
    RealType* F(int n)        { return f[n]; }

    RealType* DataF(int n)    { return dataF[n]; }
    int*      DataI(int n)    { return dataI[n]; }

    // --- Particle data

    // Number of particles
    int number;

    // Number of "ghost" particles
    int numberG;

    // Number of types of particles: 0, 1, ... , [ntypes]
    int ntypes;

    // All particles have position, force - this data is used for integration, calculating verlet lists, etc.
    RealType **x, **v, **f;
    // All particles also have a characteristic length (sg - for sigma), (inverse) mass
    RealType *sg, *im;
    // Store the type of each atom
    int *type;
    // More data, for more complex objects
    RealType **dataF; // Floating point data
    int      **dataI; // Integer data
    int       *body;  // Which body a particle is part of. No body is -1. If body == nullptr, there are no bodies

    // The actual amount of data in the arrays
    int size;

    // GFlow is a friend class -- this allows it to access protected Base members
    friend class GFlow;
  };


}
#endif // __SIM_DATA_HPP__GFLOW__