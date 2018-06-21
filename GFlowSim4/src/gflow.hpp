#ifndef __GFLOW_HPP__GFLOW__
#define __GFLOW_HPP__GFLOW__


#include "base.hpp"
#include "utility.hpp"
#include "bounds.hpp"

namespace GFlowSimulation {

  class GFlow {
  public:

    // Constructor
    GFlow(int, char**);

    // Destructor
    ~GFlow();

    // Initialize
    void initialize();

    // Run the simulation for some amount of time
    void run(RealType);

    // Data - public so anyone can access it
    class SimData *simData;             // Particle data
    class Sectorization *sectorization; // Sectorization of particles
    class Neighbors *neighbors;         // Neighbor lists
    class Communicator *communicator;   // Inter-process communicator
    class Integrator *integrator;       // Integrator

    // A vector of objects that should modify the simulation at some point(s) during execution
    vector<class Modifier*> modifiers;

  protected:
    // If true, the simulation should continue to run
    bool running;

    // How much of the requested time has been run
    RealType fulfilled_time; 

    // The number of iterations
    int iter;

    // The simulation bounds
    Bounds bounds;

    // Periodicity
    bool wrap[DIMENSIONS];

  };

}
#endif // __GFLOW_HPP__GFLOW__