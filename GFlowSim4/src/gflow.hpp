#ifndef __GFLOW_HPP__GFLOW__
#define __GFLOW_HPP__GFLOW__

#include "base.hpp"
#include "utility.hpp"

namespace GFlowSimulation {

  /*
  *  @class GFlow
  *
  *  The simulation class
  *
  */
  class GFlow {
  public:
    // Constructor
    GFlow();

    // Destructor
    ~GFlow();

    // Initialize - returns true if all pointers are non-null
    bool initialize();

    // Initialize a base object to point at GFlow's data
    void initializeBase(Base *);

    // Run the simulation for some amount of time
    void run(RealType);

    // Write data from data master to file
    void writeData(string);

    // --- Accessors
    
    // Get fulfilled time
    RealType getElapsedTime();

    // --- Mutators

    // Set all wrap values to something
    void setAllWrap(bool);

    // Creators are a friend classes --- all must be since friendship is not inherited
    friend class BoxCreator;

  protected:
    // Private helper functions
    inline void wrapPositions();

    // Data - public so anyone can access it
    class SimData *simData;             // Particle data
    class Integrator *integrator;       // Integrator
    class Sectorization *sectorization; // Sectorization of particles
    class Communicator *communicator;   // Inter-process communicator
    class DataMaster *dataMaster;       // DataMaster object for unified data collection    

    // A vector of objects that should modify the simulation at some point(s) during execution
    vector<class Modifier*> modifiers;

    // All the forces that can happen
    vector<class Force*> forces;

    // If true, the simulation should continue to run
    bool running;

    // How much of the requested time has been run
    RealType elapsed_time; 

    // The number of iterations
    int iter;

    // The simulation bounds
    Bounds bounds;

    // Periodicity
    bool wrap[DIMENSIONS];
  };

}
#endif // __GFLOW_HPP__GFLOW__