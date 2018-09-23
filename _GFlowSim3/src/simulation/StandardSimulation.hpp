/*
 * Author: Nathaniel Rupprecht
 * Start Data: January 27, 2018
 *
 */

#ifndef __STANDARD_SIMULATION_HPP__
#define __STANDARD_SIMULATION_HPP__

#include "SimulationBase.hpp"

namespace GFlow {
  
  /*
   * @class StandardSimulation
   *
   * The usual simulation class, it just runs the integrator. No fancy detections 
   * or checking
   *
   */
  class StandardSimulation : public SimulationBase {
  public:
    // Constructor
    StandardSimulation();

    // Set up the simulation
    virtual void setUp(int, char**);

    // Set parameters from the command line
    virtual void parse();

    // Run the simulation
    virtual void run();

    // Save data from the simulation
    virtual void write();    

  private:
    // Run time
    RealType runTime;
  };
  
}

#endif // __STANDARD_SIMULATION_HPP__
