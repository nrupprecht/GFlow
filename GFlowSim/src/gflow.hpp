#ifndef __GFLOW_HPP__GFLOW__
#define __GFLOW_HPP__GFLOW__

#include "base/base.hpp"
#include "utility/utility.hpp"

#include "utility/timer.hpp"

// Note to self: Using debuggers in MPI parallel: https://www.open-mpi.org/faq/?category=debugging
// mpirun -np 2 xterm -e lldb ./bin/driver

namespace GFlowSimulation {

  //! \brief What mode of running is going on.
  //!
  //! IDLE means the simulation is not currently running.
  //! INIT is the mode for initialization, like relaxing objects.
  //! SIM is the mode for actually running a simulation.
  enum class RunMode { IDLE, INIT, SIM };

  /*
  *  \brief GFlow
  *
  *  The simulation class.
  *
  */
  class GFlow {
  public:
    //! \brief Constructor.
    GFlow(int);

    //! \brief Destructor.
    ~GFlow();

    //! \brief Initialize - returns true if all pointers are non-null.
    bool initialize();

    //! \brief Initialize a base object to point at GFlow's data.
    void initializeBase(Base*);

    //! \brief Run the simulation for some amount of time.
    void run(long double=-1);

    //! \brief Write data from data master to file.
    void writeData(string);

    // --- Accessors

    //! \brief Whether a SimData has been allocated.
    bool hasSimData()     { return simData!=nullptr; }
    //! \brief Whether an integrator has been allocated.
    bool hasIntegrator()  { return integrator!=nullptr; }
    //! \brief Whether an interaction handler has been allocated.
    bool hasHandler()      { return handler!=nullptr; }
    //! \brief Whether a data master has been allocated.
    bool hasDataMaster()  { return dataMaster!=nullptr; }
    //! \brief Whether a force master has been allocated.
    bool hasForceMaster() { return forceMaster!=nullptr; }

    //! \brief Get argc
    int getArgC() const { return argc; }

    //! \brief Get argv
    char** getArgV() const { return argv; }

    //! \brief Get the requested time.
    long double getRequestedTime() const;

    //! \brief Get the total amount of time ever requested.
    //!
    //! Get the total amount of time that has ever been requested of this
    //! GFlow object. Time can be requested and run, and then more time can
    //! be requested and run, so the most recent amount of requested time is
    //! not necessarily the total amount of time ever requested.
    long double getTotalRequestedTime() const;
    
    //! \brief Get fulfilled time.
    long double getElapsedTime() const;

    //! \brief Get all the time the simulation ran for.
    long double getTotalTime() const;

    //! \brief Get the strength of the boundary force.
    RealType getBoundaryForce() const;

    //! \brief Get the energy associated with the domain.
    RealType getBoundaryEnergy() const;

    //! \brief Get the current time step.
    RealType getDT() const;

    //! \brief Get the iteration.
    long int getIter() const;

    //! \brief Get the number of forces
    int getNumInteractions() const;

    //! \brief Get the number of particles
    int getNumParticles() const;

    //! \brief Get the bounds
    Bounds getBounds() const;

    //! \brief Get the boundary conditions
    const BCFlag* getBCs() const;

    //! \brief Get a single boundary condition
    BCFlag getBC(int) const;

    //! \brief Get the number of types of particles in the simulation
    int getNTypes() const;

    //! \brief Get the dimensionality of the simulation.
    int getSimDimensions() const;

    //! \brief Get the use forces flag.
    bool getUseForces() const;

    //! \brief Get the command line arguments.
    pair<int, char**> getCommand() const;

    //! \brief Get the vector of nonbonded interactions.
    const vector<class Interaction*>& getInteractions() const;

    //! \brief Get the vector of bonded interactions.
    const vector<class Bonded*>& getBondedInteractions() const;

    //! \brief Get the sim data object.
    class SimData* getSimData();

    //! \brief Get the data master object.
    class DataMaster* getDataMaster();

    //! \brief Get the force master object.
    class ForceMaster* getForceMaster() ;

    //! \brief Get the integrator.
    class Integrator* getIntegrator();

    //! \brief Get the number of integrators that gflow has.
    int getNumIntegrators() const; 

    //! \brief Get the minimum image displacement between two positions.
    void getDisplacement(const RealType*, const RealType*, RealType*);

    //! \brief Get the minimum image of a displacement vector.
    void minimumImage(RealType*);

    //! \brief Get the minimum image distance of a single component.
    void minimumImage(RealType&, int);

    RealType getDistance(const RealType*, const RealType*);

    //! \brief Get the run mode of the simulation.
    RunMode getRunMode();

    // --- Mutators

    //! \brief Add an interaction. 
    //!
    //! GFlow only adds the interaction if it is non null.
    void addInteraction(Interaction*);

    //! \brief Add a bonded interaction.
    void addBonded(class Bonded*);

    //! \brief Add a body.
    void addBody(class Body*);

    //! \brief Add another integrator.
    void addIntegrator(class Integrator*);

    //! \brief Set the command info
    void setCommand(int, char**);

    //! \brief Set all wrap values to the same value.
    void setAllBCs(BCFlag);

    //! \brief Set a single boundary condition.
    void setBC(const int, const BCFlag);

    //! \brief Set the use forces flag.
    void setUseForces(bool);

    //! \brief Set the simulation bounds.
    void setBounds(const Bounds&);

    //! \brief Set the repulsion stength for repulsing boundary conditions.
    void setRepulsion(RealType);

    //! \brief Set the dissipation stength for repulsing boundary conditions.
    void setDissipation(RealType);

    //! \brief Set the attraction acceleration.
    void setAttraction(RealType);

    //! \brief Set the running flag.
    void setRunning(bool);

    //! \brief Set the print updates flag.
    void setPrintUpdates(bool);

    //! \brief Set the update printing interval.
    void setUpdateInterval(RealType);

    //! \brief Set the run mode.
    void setRunMode(RunMode);

    //! \brief Set the amount of time we should run for.
    void requestTime(RealType);

    //! \brief Set the elapsed time.
    void setElapsedTime(RealType);

    //! \brief Keep positions in bounds.
    void wrapPositions();

    //! \brief "Reflect" positions off the bounds.
    void reflectPositions();

    //! \brief Apply a force to keep particles in bounds.
    void repulsePositions();

    //! \brief Apply a harmonic force to keep particles attracted to the center of the simulation
    void attractPositions();

    //! \brief Instructs the interaction handler to remove particles that are overlapping by more than some fraction.
    void removeOverlapping(RealType);

    //! \brief Add a data object.
    void addDataObject(class DataObject*);

    //! \brief Add a modifier object.
    void addModifier(class Modifier*);

    //! \brief Reset all timers (use e.g. after doing relaxation of a random initial state).
    void resetAllTimes();

    //! \brief Set the start recording time.
    void setStartRecTime(RealType);

    //! \brief Set the frames per second for all data objects.
    void setFPS(RealType);

    //! \brief Set the fps of particular data objects.
    void setFPS(int, RealType);

    //! \brief Set the time step.
    void setDT(RealType);

    //! \brief Set the max timestep.
    void setMaxDT(RealType);

    //! \brief Set data master command line data.
    void setDMCmd(int, char**);

    void giveFileToDataMaster(string, string);

    //! \brief Get the boltzmann constant.
    RealType getKB();

    // Creators are a friend classes --- all must be since friendship is not inherited
    friend class Creator;
    friend class BoxCreator;
    friend class BondBoxCreator;
    friend class BipartiteBoxCreator;
    friend class DebugCreator;
    friend class FlowCreator;
    friend class FileParseCreator;
    friend class LineCrystalCreator;
    friend class PolymerCreator;
    friend class FillAreaCreator;

    // Other friend classes
    friend class ForceMaster;
    friend class DataMaster;
    friend class DataObject;
    friend class Base;

  protected:
    // --- Private helper functions
    //! Clear all the force arrays.
    inline void clearForces();

    //! \brief Check whether any modifiers need to be removed.
    inline void handleModifiers();

    // --- Data - public so anyone can access it
    class SimData     *simData = nullptr;      // Particle data
    class Integrator  *integrator = nullptr;   // Integrator
    class InteractionHandler *handler = nullptr;      
    class DataMaster  *dataMaster = nullptr;   // DataMaster object for unified data collection  
    class ForceMaster *forceMaster = nullptr;  // ForceMaster object for defining and storing interparticle forces  
    class Topology    *topology = nullptr;     // Processor topology

    //! \brief Additional integrators. \todo Do this better.
    vector<class Integrator*> additional_integrators;

    //! \brief A vector of objects that should modify the simulation at some point(s) during execution.
    std::list<class Modifier*> modifiers;

    //! \brief All the short range, non-bonded, forces that can happen - which ones correspond to which pairs of particles is controlled by
    // the ForceMaster object.
    vector<class Interaction*> interactions;

    //! \brief All the bonded forces that can happen.
    vector<class Bonded*> bondedInteractions;

    //! \brief All the bodies in the simulation.
    vector<class Body*> bodies;

    //! \brief If true, the simulation should continue to run.
    bool running;

    //! \brief If true, do tasks related to force computation.
    bool useForces = true;

    //! \brief How much time we have been requested to run for.
    long double requested_time;

    //! \brief How much time has ever been requested.
    long double total_requested_time;

    //! \brief How much of the requested time has been run.
    long double elapsed_time; 

    //! \brief How much time has been run over all runs.
    long double total_time;

    //! \brief The number of iterations that have passed.
    long int iter;

    //! \brief The simulation bounds.
    Bounds bounds;

    //! \brief Boundary types.
    BCFlag *boundaryConditions;

    //! \brief The number of dimensions
    int sim_dimensions;

    //! \brief Strength of boundary repulsion forces.
    RealType repulsion;

    //! \brief Dissipation for the boundary repulsion forces.
    RealType dissipation;

    //! \brief The attraction towards the center of the simulation
    RealType center_attraction;

    //! \brief Total boundary force applied this iteration.
    RealType boundaryForce;

    //! \brief Energy due to e.g. particles being repulsed by a boundary potential.
    RealType boundaryEnergy = 0;

    // The command info (optional)
    int argc;
    char **argv;

    // Timing
    Timer domain_timer;

    //! \brief Timer to time how long the bonded interactions take.
    Timer bonded_timer;

    //! \brief Timer to time how long the body execution takes.
    Timer body_timer;

    //! \brief Whether to print updates to a screen.
    bool print_updates = false;

    //! \brief The ostream with which to print updates.
    std::ostream *monitor = &cout;

    //! \brief At what intervals to print updates (in simulation seconds).
    RealType update_interval = 250.;

    //! \brief The run mode of the simulation.
    RunMode runMode = RunMode::IDLE;

    // --- Physical constants. These are static, since they should be the same for all objects

    //! \brief Boltzmann constant.
    RealType KB = 1.; 
  };

}
#endif // __GFLOW_HPP__GFLOW__
