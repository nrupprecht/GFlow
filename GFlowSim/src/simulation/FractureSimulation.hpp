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
    inline void checkForBreaks();
    inline void compare(const VListSubType&, const VListSubType&, vector<pair<int,int> >&, vector<pair<int,int> >&);

    // Run time
    RealType runTime;

    // A record of the verlet lists
    VListType verletListRecord;
    
    // A lists {time, particle #, particle #} of when breakages and bondings occur
    vector<std::tuple<RealType,int,int> > breakages, bondings;
    // Cumulative number of breaks and bonds as a time series
    vector<pair<RealType,int> > cumNumBreaks, cumNumBonds;
    // Cumulative number of breaks and bonds
    int numBreaks, numBonds;
  };

}

#endif // __FRACTURE_SIMULATION_HPP__
