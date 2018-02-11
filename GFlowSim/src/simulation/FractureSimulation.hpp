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

    // Get breakage, bonding data
    auto& getBreakages() { return breakages; }
    auto& getBondings()  { return bondings; }
    auto& getCumNumBreaks() { return cumNumBreaks; }
    auto& getCumNumBonds()  { return cumNumBonds; }

  private:
    // Private helper functions
    inline void checkForBreaks();
    inline void compare(vector<pair<int,int> >&, vector<pair<int,int> >&);
    inline void setVerletListRecord(VListType&);
    inline void setVerletListRecord();
    inline void strengthenProcedure();
    inline RealType getHardCoreDistance(int, int);

    // Run time
    RealType limitTime;

    // Actual time simulation took
    RealType totalTime;

    // Delay between checking for breaks/bonds
    RealType delay;

    // A sectorization for doing breakage analysis
    Sectorization bondLists;

    // A record of the verlet lists for the large neighborhood
    VListType verletListRecord;

    // An unordered set of all the particles that have most recently broken bonds
    std::unordered_set<int> marked;

    // A way to sort the verlet list we get from the sectorization
    vector<VListSubType*> accessor;

    // Particles that were closer than the inner radius
    vector<vector<int> > bondedParticles; 

    // Radius required for bonding - distance the rmin points should be from each other
    RealType bondRadius;

    // Whether to mark breakages or the formation of bonds
    bool markBreaks;

    // Time that last break occured
    RealType lastBreak;

    // Amount of time we should animate breaks
    RealType breakDisplayLength;

    // Whether to do iterative strengthening
    bool strengthen;

    // A list of [trial][break time]
    vector<vec2> breakTimes;

    // Slot in dataRecord that we use to record breakTimes and improvement
    int breakDataSlot, improvementSlot;

    // The strengthening iteration we are on and the total number of strengthening iterations
    int heatIters, maxIters;
    
    // How long to apply heat and relax during a strengthening
    RealType heatTime, relaxTime;

    // The temperature to use for heating
    RealType heatTemperature;

    // The radius of the heat application
    RealType heatRadius;
    
    // The dF/dt of the wall
    RealType dFx, dFy;

    // A lists {time, particle #, particle #} of when breakages and bondings occur
    vector<std::tuple<RealType,int,int> > breakages, bondings;
    // Cumulative number of breaks and bonds as a time series
    vector<pair<RealType,int> > cumNumBreaks, cumNumBonds;
    // Cumulative number of breaks and bonds
    int numBreaks, numBonds;

    // Force to cause the first break
    RealType breakF0;
  };

}

#endif // __FRACTURE_SIMULATION_HPP__
