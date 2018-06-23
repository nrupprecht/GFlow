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
    void run(RealType=-1);

    // Write data from data master to file
    void writeData(string);

    // --- Accessors
    
    // Get fulfilled time
    RealType getElapsedTime() const;

    // Get the bounds
    Bounds getBounds() const;

    // Get the wrapping data
    const bool* getWrap() const;

    // --- Mutators

    // Set all wrap values to something
    void setAllWrap(bool);

    // Set the amount of time we should run for
    void requestTime(RealType);

    // Keep positions in bounds
    void wrapPositions();

    // Creators are a friend classes --- all must be since friendship is not inherited
    friend class BoxCreator;

    // Force master is a friend class
    friend class ForceMaster;

  protected:
    // Private helper functions

    // Data - public so anyone can access it
    class SimData *simData;             // Particle data
    class Integrator *integrator;       // Integrator
    class Sectorization *sectorization; // Sectorization of particles
    class Communicator *communicator;   // Inter-process communicator
    class DataMaster *dataMaster;       // DataMaster object for unified data collection  
    class ForceMaster *forceMaster;     // ForceMaster object for defining and storing interparticle forces  

    // A vector of objects that should modify the simulation at some point(s) during execution
    vector<class Modifier*> modifiers;

    // All the forces that can happen
    vector<class Force*> forces;

    // If true, the simulation should continue to run
    bool running;

    // How much time we have been requested to run for
    RealType requested_time;

    // How much of the requested time has been run
    RealType elapsed_time; 

    // How much time has been run over all runs
    RealType total_time;

    // The number of iterations
    int iter;

    // The simulation bounds
    Bounds bounds;

    // Periodicity
    bool wrap[DIMENSIONS];
  };

}
#endif // __GFLOW_HPP__GFLOW__
