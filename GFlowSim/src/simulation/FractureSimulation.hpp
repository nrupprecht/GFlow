/*
 * Author: Nathaniel Rupprecht
 * Start Data: January 27, 2018
 *
 */

#ifndef __FRACTURE_SIMULATION_HPP__
#define __FRACTURE_SIMULATION_HPP__

#include "SimulationBase.hpp"

namespace GFlow {

  class FractureSimulation : public SimulationBase {
  public:
    // Constructor
    FractureSimulation();

    virtual void setUp(int, char**);

    // Set parameters from the command line
    virtual void parse();

    // Run the simulation
    virtual void run();

    // Save data from the simulation
    virtual void write();

  private:
    // Private helper functions
    void checkForBreaks();

    // Run time
    RealType runTime;
  };

}

#endif // __FRACTURE_SIMULATION_HPP__
