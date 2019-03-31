#ifndef __GFLOW_HPP__GFLOW__
#define __GFLOW_HPP__GFLOW__

#include "base/base.hpp"
#include "utility/utility.hpp"

#include "utility/timer.hpp"

namespace GFlowSimulation {

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
    void initializeBase(Base *);

    //! \brief Run the simulation for some amount of time.
    void run(RealType=-1);

    //! \brief Write data from data master to file.
    void writeData(string);

    // --- Accessors

    //! \brief Whether a SimData has been allocated.
    bool hasSimData()     { return simData!=nullptr; }
    //! \brief Whether an integrator has been allocated.
    bool hasIntegrator()  { return integrator!=nullptr; }
    //! \brief Whether a domain has been allocated.
    bool hasDomain()      { return domain!=nullptr; }
    //! \brief Whether a data master has been allocated.
    bool hasDataMaster()  { return dataMaster!=nullptr; }
    //! \brief Whether a force master has been allocated.
    bool hasForceMaster() { return forceMaster!=nullptr; }

    //! \brief Get argc
    int getArgC() const { return argc; }

    //! \brief Get argv
    char** getArgV() const { return argv; }

    //! \brief Get the requested time.
    RealType getRequestedTime() const;

    //! \brief Get the total amount of time ever requested.
    //!
    //! Get the total amount of time that has ever been requested of this
    //! GFlow object. Time can be requested and run, and then more time can
    //! be requested and run, so the most recent amount of requested time is
    //! not necessarily the total amount of time ever requested.
    RealType getTotalRequestedTime() const;
    
    //! \brief Get fulfilled time.
    RealType getElapsedTime() const;

    //! \brief Get all the time the simulation ran for.
    RealType getTotalTime() const;

    //! \brief Get the strength of the boundary force.
    RealType getBoundaryForce() const;

    //! \brief Get the current time step.
    RealType getDT() const;

    //! \brief Get the iteration.
    int getIter() const;

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

    int getSimDimensions() const;

    pair<int, char**> getCommand() const;

    const vector<class Interaction*>& getInteractions() const;

    class SimData* getSimData();

    class DataMaster* getDataMaster();

    class ForceMaster* getForceMaster() ;

    class Integrator* getIntegrator();

    const RealType* getVComCorrection() const;

    void getDisplacement(const RealType*, const RealType*, RealType*);

    RealType getDistance(const RealType*, const RealType*);

    // --- Mutators

    //! \brief Add an interaction. 
    //!
    //! GFlow only adds the interaction if it is non null.
    void addInteraction(Interaction*);

    //! \brief Set the command info
    void setCommand(int, char**);

    //! \brief Set the bounds
    void setBounds(Bounds);

    //! \brief Set all wrap values to the same value.
    void setAllBCs(BCFlag);

    //! \brief Set a single boundary condition.
    void setBC(const int, const BCFlag);

    //! \brief Set the repulsion stength for repulsing boundary conditions.
    void setRepulsion(RealType);

    //! \brief Set the dissipation stength for repulsing boundary conditions.
    void setDissipation(RealType);

    //! \brief Set the attraction acceleration.
    void setAttraction(RealType);

    //! \brief Set the running flag.
    void setRunning(bool);

    //! \brief Set the amount of time we should run for.
    void requestTime(RealType);

    //! \brief Keep positions in bounds.
    void wrapPositions();

    //! \brief "Reflect" positions off the bounds.
    void reflectPositions();

    //! \brief Apply a force to keep particles in bounds.
    void repulsePositions();

    //! \brief Apply a harmonic force to keep particles attracted to the center of the simulation
    void attractPositions();

    //! \brief Instructs the domain to remove particles that are overlapping by more than some fraction.
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

    //! \brief Set data master command line data.
    void setDMCmd(int, char**);

    void giveFileToDataMaster(string, string);

    // Creators are a friend classes --- all must be since friendship is not inherited
    friend class Creator;
    friend class BoxCreator;
    friend class BondBoxCreator;
    friend class BipartiteBoxCreator;
    friend class DebugCreator;
    friend class FlowCreator;
    friend class FileParseCreator;
    friend class LineCrystalCreator;

    // Force master is a friend class
    friend class ForceMaster;

    // Data master is a friend class
    friend class DataMaster;

    friend class Base;

    friend class PolymerCreator;

  protected:
    // --- Private helper functions
    //! Clear all the force arrays.
    inline void clearForces();

    //! \brief Check whether any modifiers need to be removed.
    inline void handleModifiers();

    // --- Data - public so anyone can access it
    class SimData     *simData = nullptr;      // Particle data
    class Integrator  *integrator = nullptr;   // Integrator
    class DomainBase  *domain = nullptr;       // Domain
    class DataMaster  *dataMaster = nullptr;   // DataMaster object for unified data collection  
    class ForceMaster *forceMaster = nullptr;  // ForceMaster object for defining and storing interparticle forces  
    class Topology    *topology = nullptr;     // Processor topology

    //! \brief A vector of objects that should modify the simulation at some point(s) during execution.
    std::list<class Modifier*> modifiers;

    //! \brief All the forces that can happen - which ones correspond to which pairs of particles is controlled by
    // the ForceMaster object.
    vector<class Interaction*> interactions;

    //! \brief If true, the simulation should continue to run.
    bool running;

    //! \brief If true, do tasks related to force computation.
    bool useForces;

    //! \brief How much time we have been requested to run for.
    RealType requested_time;

    //! \brief How much time has ever been requested.
    RealType total_requested_time;

    //! \brief How much of the requested time has been run.
    RealType elapsed_time; 

    //! \brief How much time has been run over all runs.
    RealType total_time;

    //! \brief The number of iterations that have passed.
    int iter;

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

    //! Total boundary force applied this iteration.
    RealType boundaryForce;

    // The command info (optional)
    int argc;
    char **argv;

    // Timing
    Timer fhs_timer;
    Timer shs_timer;
    Timer domain_timer;
    Timer forces_timer;

    // MPI
    int rank = 0;
    int numProc = 1;
  };

}
#endif // __GFLOW_HPP__GFLOW__
