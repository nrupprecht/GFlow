/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __INTEGRATOR_HPP__
#define __INTEGRATOR_HPP__

// Includes
#include "../control/SimData.hpp"
#include "../data/DataRecord.hpp"
 
namespace GFlow {

  /*
   * @class SimData
   * Base class for integrators which use forces/torques to advance the simulation
   */
  class Integrator {
  public:
    // Constructor
    Integrator();

    // Initialization
    void initialize();

    // Integrate function wraps protected virtual function
    void integrate() {
      _integrate();
    }

  protected:
    // Private virtual functions
    virtual void _integrate() {} ; // Should be purely abstract

    // Time step
    double dt;
    // Current time
    double time;
    // Simulation time
    double runTime;

    // Current iteration
    int iter;
    // Max iteration
    int maxIter;
    
    // Pointer to the simulation data
    SimData* simData;

    // Pointer to data recorder
    DataRecord* dataRecord;
    
  };

}
#endif // __INTEGRATOR_HPP__
