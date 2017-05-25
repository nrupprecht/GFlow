/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __INTEGRATOR_HPP__
#define __INTEGRATOR_HPP__

// Includes
#include "../control/SimData.hpp"
#include "../control/Sectorization.hpp"
#include "../data/DataRecord.hpp"
#include "../control/ForceHandler.hpp"
#include "../control/StatusObject.hpp"
#include "../control/TerminationCondition.hpp"
 
namespace GFlow {

  // Forward declaration to DataRecord
  class DataRecord;

  /*
   * @class SimData
   * Base class for integrators which use forces/torques to advance the simulation
   */
  class Integrator {
  public:
    // Default constructor
    Integrator();

    // SimData initializing constructor
    Integrator(SimData*);

    // SimData and DataRecord initializing constructor
    Integrator(SimData*, DataRecord*);

    // Destructor - cleans up sectors
    ~Integrator();

    // Initialization -- Add the run time here for now
    void initialize(RealType runTime);

    // Give the integrator a data record object
    void setDataRecord(DataRecord* dr);

    // Integrate function wraps protected virtual function
    void integrate();

    // Get the timestep
    RealType getDt() { return dt; }

    // Set the timestep
    void setDt(RealType t) { dt = t; }

    // Add a termination condition - Integrator is responsible for managing them
    void addTerminationCondition(TerminationCondition *t) { termination.push_back(t); }

    // DataRecord is a friend class
    friend class DataRecord;

  protected:
    // Private virtual functions
    virtual void initializeSectors();
    virtual void initializeForceHandler();
    virtual void preIntegrate();
    virtual inline void _integrate() = 0;
    virtual void postIntegrate();
    virtual void preStep();
    virtual void postStep();

    // Whether the simulation is running or not
    bool running;

    // Time step
    RealType dt;

    // Current time
    RealType time;

    // How much simulated time the simulation should run for (barring some early termination condition)
    RealType runTime;

    // Current iteration
    int iter;
    
    // Pointer to the simulation data
    SimData* simData;

    // Sectorization object
    Sectorization *sectors;

    // Pointer to data recorder
    DataRecord* dataRecord;

    // Force handler
    ForceHandler* forceHandler;

    // List of terminal conditions - Integrator is responsible for managing them
    vector<TerminationCondition*> termination;
  };

}
#endif // __INTEGRATOR_HPP__
