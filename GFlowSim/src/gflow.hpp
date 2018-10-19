#ifndef __GFLOW_HPP__GFLOW__
#define __GFLOW_HPP__GFLOW__

#include "base/base.hpp"
#include "utility/utility.hpp"

namespace GFlowSimulation {

  /*
  *  @brief GFlow
  *
  *  The simulation class
  *
  */
  class GFlow {
  public:
    //! Constructor.
    GFlow();

    //! Destructor.
    ~GFlow();

    //! Initialize - returns true if all pointers are non-null.
    bool initialize();

    //! Initialize a base object to point at GFlow's data.
    void initializeBase(Base *);

    //! Run the simulation for some amount of time.
    void run(RealType=-1);

    //! Write data from data master to file.
    void writeData(string);

    // --- Accessors

    //! @brief Whether a SimData has been allocated.
    bool hasSimData()     { return simData!=nullptr; }
    //! @brief Whether an integrator has been allocated.
    bool hasIntegrator()  { return integrator!=nullptr; }
    //! @brief Whether a domain has been allocated.
    bool hasDomain()      { return domain!=nullptr; }
    //! @brief Whether a data master has been allocated.
    bool hasDataMaster()  { return dataMaster!=nullptr; }
    //! @brief Whether a force master has been allocated.
    bool hasForceMaster() { return forceMaster!=nullptr; }

    //! Get argc
    int getArgC() const { return argc; }

    //! @brief Get argv
    char** getArgV() const { return argv; }

    //! Get the requested time.
    RealType getRequestedTime() const;

    //! @brief Get the total amount of time ever requested.
    //!
    //! Get the total amount of time that has ever been requested of this
    //! GFlow object. Time can be requested and run, and then more time can
    //! be requested and run, so the most recent amount of requested time is
    //! not necessarily the total amount of time ever requested.
    RealType getTotalRequestedTime() const;
    
    //! Get fulfilled time.
    RealType getElapsedTime() const;

    //! Get all the time the simulation ran for.
    RealType getTotalTime() const;

    //! Get the strength of the boundary force.
    RealType getBoundaryForce() const;

    //! Get the current time step.
    RealType getDT() const;

    //! Get the number of iterations the simulation had run for.
    int getIter() const;

    // Get the number of forces
    int getNumInteractions() const;

    // Get the number of particles
    int getNumParticles() const;

    // Get the bounds
    Bounds getBounds() const;

    // Get the boundary conditions
    const BCFlag* getBCs() const;

    // Get a single boundary condition
    BCFlag getBC(int) const;

    // Get the number of types of particles in the simulation
    int getNTypes() const;

    pair<int, char**> getCommand() const;

    const vector<class Interaction*>& getInteractions() const;

    class DataMaster* getDataMaster() const;

    // --- Mutators

    //! Set the command info
    void setCommand(int, char**);

    //! Set the bounds
    void setBounds(Bounds);

    //! Set all wrap values to the same value.
    void setAllBCs(BCFlag);

    //! @brief Set a single boundary condition.
    void setBC(const int, const BCFlag);

    //! Set the repulsion stength for repulsing boundary conditions.
    void setRepulsion(RealType);

    //! Set the dissipation stength for repulsing boundary conditions.
    void setDissipation(RealType);

    //! Set the amount of time we should run for.
    void requestTime(RealType);

    //! Keep positions in bounds.
    void wrapPositions();

    //! "Reflect" positions off the bounds.
    void reflectPositions();

    //! Apply a force to keep particles in bounds.
    void repulsePositions();

    //! Add a data object.
    void addDataObject(class DataObject*);

    //! Add a modifier object.
    void addModifier(class Modifier*);

    //! Reset all timers (use e.g. after doing relaxation of a random initial state).
    void resetAllTimes();

    //! Set the start recording time.
    void setStartRecTime(RealType);

    //! Set the frames per second for all data objects.
    void setFPS(RealType);

    //! Set the fps of particular data objects.
    void setFPS(int, RealType);

    //! Set the time step.
    void setDT(RealType);

    //! Set data master command line data.
    void setDMCmd(int, char**);

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

  protected:
    // --- Private helper functions
    //! Clear all the force arrays.
    inline void clearForces();

    //! @brief Check whether any modifiers need to be removed.
    inline void handleModifiers();

    // --- Data - public so anyone can access it
    class SimData *simData;             // Particle data
    class Integrator *integrator;       // Integrator
    class DomainBase *domain;           // Domain
    class DataMaster *dataMaster;       // DataMaster object for unified data collection  
    class ForceMaster *forceMaster;     // ForceMaster object for defining and storing interparticle forces  

    class BondData  *bondData;
    class AngleData *angleData;

    //! @brief A vector of objects that should modify the simulation at some point(s) during execution.
    std::list<class Modifier*> modifiers;

    //! @brief All the forces that can happen - which ones correspond to which pairs of particles is controlled by
    // the ForceMaster object.
    vector<class Interaction*> interactions;

    //! @brief If true, the simulation should continue to run.
    bool running;

    //! @brief If true, do tasks related to force computation.
    bool useForces;

    //! @brief How much time we have been requested to run for.
    RealType requested_time;

    //! @brief How much time has ever been requested.
    RealType total_requested_time;

    //! @brief How much of the requested time has been run.
    RealType elapsed_time; 

    //! @brief How much time has been run over all runs.
    RealType total_time;

    //! @brief The number of iterations that have passed.
    int iter;

    //! @brief The simulation bounds.
    Bounds bounds;

    //! @brief Boundary types.
    BCFlag boundaryConditions[DIMENSIONS];

    //! @brief The number of dimensions
    int sim_dimensions;

    //! @brief Strength of boundary repulsion forces.
    RealType repulsion;

    //! @brief Dissipation for the boundary repulsion forces.
    RealType dissipation;

    //! Total boundary force applied this iteration.
    RealType boundaryForce;

    // The command info (optional)
    int argc;
    char **argv;
  };

}
#endif // __GFLOW_HPP__GFLOW__
