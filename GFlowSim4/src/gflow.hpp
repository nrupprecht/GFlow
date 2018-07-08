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

    RealType getBoundaryForce() const;

    RealType getDT() const;

    // Get the number of iterations
    int getIter() const;

    // Get the number of forces
    int getNumForces() const;

    // Get the bounds
    Bounds getBounds() const;

    // Get the boundary conditions
    const BCFlag* getBCs() const;

    // Get a single boundary condition
    BCFlag getBC(int) const;

    // Get the number of types of particles in the simulation
    int getNTypes() const;

    pair<int, char**> getCommand() const;

    const vector<class Force*>& getForces() const;

    // --- Mutators

    // Set the command info
    void setCommand(int, char**);

    // Set the bounds
    void setBounds(Bounds);

    // Set all wrap values to something
    void setAllBCs(BCFlag);

    void setRepulsion(RealType);

    // Set the amount of time we should run for
    void requestTime(RealType);

    // Keep positions in bounds
    void wrapPositions();

    // "Reflect" positions off the bounds
    void reflectPositions();

    // Apply a force to keep particles in bounds
    void repulsePositions();

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

    // Set the time step
    void setDT(RealType);

    // Set data master command line data
    void setDMCmd(int, char**);

    // Creators are a friend classes --- all must be since friendship is not inherited
    friend class Creator;
    friend class BoxCreator;
    friend class BondBoxCreator;
    friend class BipartiteBoxCreator;
    friend class DebugCreator;

    // Force master is a friend class
    friend class ForceMaster;

    // Data master is a friend class
    friend class DataMaster;

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

    // Boundary types
    BCFlag boundaryConditions[DIMENSIONS];

    // Strength of boundary repulsion forces, and total force applied
    RealType repulsion, boundaryForce;

    // The command info (optional)
    int argc;
    char **argv;
  };

}
#endif // __GFLOW_HPP__GFLOW__
