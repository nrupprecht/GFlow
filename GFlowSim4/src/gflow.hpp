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

    // Get the requested time
    RealType getRequestedTime() const;

    RealType getTotalRequestedTime() const;
    
    // Get fulfilled time
    RealType getElapsedTime() const;

    // Get all the time the simulation ran for
    RealType getTotalTime() const;

    // Get the number of iterations
    int getIter() const;

    // Get the number of forces
    int getNumForces() const;

    // Get the bounds
    Bounds getBounds() const;

    // Get the wrapping data
    const bool* getWrap() const;

    // Get the number of types of particles in the simulation
    int getNTypes() const;

    pair<int, char**> getCommand() const;

    // --- Mutators

    // Set the command info
    void setCommand(int, char**);

    // Set the bounds
    void setBounds(Bounds);

    // Set all wrap values to something
    void setAllWrap(bool);

    // Set the amount of time we should run for
    void requestTime(RealType);

    // Keep positions in bounds
    void wrapPositions();

    // Add a data object
    void addDataObject(class DataObject*);

    // Add a modifier
    void addModifier(class Modifier*);

    // Reset all timers (use e.g. after doing relaxation of a random initial state)
    void resetAllTimes();

    // Set the start recording time
    void setStartRecTime(RealType);

    // Set the frames per second for all data objects
    void setFPS(RealType);

    // Set the fps of particular data objects
    void setFPS(int, RealType);

    // Creators are a friend classes --- all must be since friendship is not inherited
    friend class Creator;
    friend class BoxCreator;
    friend class BondBoxCreator;
    friend class BinaryBoxCreator;
    friend class DebugCreator;

    // Force master is a friend class
    friend class ForceMaster;

  protected:
    // --- Private helper functions
    inline void clearForces();

    // --- Data - public so anyone can access it
    class SimData *simData;             // Particle data
    class Integrator *integrator;       // Integrator
    class Sectorization *sectorization; // Sectorization of particles
    class Communicator *communicator;   // Inter-process communicator
    class DataMaster *dataMaster;       // DataMaster object for unified data collection  
    class ForceMaster *forceMaster;     // ForceMaster object for defining and storing interparticle forces  

    // A vector of objects that should modify the simulation at some point(s) during execution
    vector<class Modifier*> modifiers;

    // All the forces that can happen - which ones correspond to which pairs of particles is controlled by
    // the ForceMaster object
    vector<class Force*> forces;

    // If true, the simulation should continue to run
    bool running;

    // How much time we have been requested to run for
    RealType requested_time;

    // How much time has ever been requested
    RealType total_requested_time;

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

    // The command info (optional)
    int argc;
    char **argv;
  };

}
#endif // __GFLOW_HPP__GFLOW__
