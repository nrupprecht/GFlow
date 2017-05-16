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
 
namespace GFlow {

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
    void initialize(double runTime);

    // Give the integrator a data record object
    void setDataRecord(DataRecord* dr);

    // Integrate function wraps protected virtual function
    void integrate();

  protected:
    // Private virtual functions, purely abstract
    virtual void _integrate() = 0;

    // Private virtual functions
    virtual void preIntegrate();
    virtual void postIntegrate();

    // Time step
    RealType dt;

    // Current time
    RealType time;

    // How much simulated time the simulation should run for
    RealType runTime;

    // Current iteration
    int iter;

    // Max iteration
    int maxIter;
    
    // Pointer to the simulation data
    SimData* simData;

    // Pointer to sectorization
    Sectorization* sectors;

    // Pointer to data recorder
    DataRecord* dataRecord;

    // Force handler
    ForceHandler* forceHandler;
    
  };

}
#endif // __INTEGRATOR_HPP__
