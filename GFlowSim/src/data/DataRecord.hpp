/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __DATA_RECORD_HPP__
#define __DATA_RECORD_HPP__

// Includes
#include "../../include/CSVUtility.hpp"
#include "../control/SimData.hpp"
#include "../control/Sectorization.hpp"
#include "StatFunc.hpp"

namespace GFlow {

  // Forward declaration to Integrator
  class Integrator;

  /*
   * @class DataRecord
   * This class generates statistics from simulation data
   *
   */
  class DataRecord {
  public:
    // Default constructor
    DataRecord();

    // Start and end timing
    void startTiming();
    void endTiming();
   
    // Set lastRecord time for performance timing
    void markTime();

    // Get the elapsed time
    double getElapsedTime() const;

    // (Potentially) record data related to simulation data
    void record(SimData*, RealType);

    // Get sectorization data at the end
    void getSectorizationData(Sectorization*);

    // Write data
    void writeData(SimData* = nullptr) const;

    // Write run summary data
    void writeRunSummary(SimData* = nullptr, Integrator* = nullptr) const;

    /*** Mutators ***/

    // Set options
    void setRecPos(bool b)  { recPos = b; }
    void setRecPerf(bool b) { recPerf = b; }

    // Set the write directory
    void setWriteDirectory(string w) { writeDirectory = w; }

    // Set the (expected) run time
    void setRunTime(RealType t) { runTime = t; }

    // Set the actual simulated time
    void setActualTime(RealType t) { actualTime = t; }

    // Add a stat function
    void addStatFunction(StatFunc, string);

    /*** Accessors ***/

    // Get the simulation run time
    RealType getRunTime() const { return runTime; }

    // Get the number of stat functions
    int getNumberOfStatFunctions() const { return statFunctionData.size(); }

    // Get the data recorded by the stat function
    vector<pair<RealType, RealType> > getStatFunctionData(int) const;
    
    // Get the name of the stat function
    string getStatFunctionName(int) const;

    // Get the performance record
    vector<pair<RealType, RealType> > getPerformanceRecord() const;

    /*** Exception classes ***/
    struct BadStatFunction {
      BadStatFunction(int i) : value(i) {};
      int value;
    };

  private:
    // Private helper functions
    void writeParticleData(std::ofstream&, SimData*) const;

    // The name of the directory we should write to 
    string writeDirectory;

    // How long between recording data
    RealType delay;

    // What time we last recorded data
    RealType lastRecord;

    // How many iterations we have recorded
    int recIter;

    // How long the simulation should last and how long the simulation actually lasted
    RealType runTime, actualTime;

    // Sectorization information
    int nsx, nsy;
    RealType sdx, sdy;
    RealType cutoff, skinDepth;
    int numberOfVerletLists;
    RealType avePerVerletList;
    int occupiedSectors;
    RealType avePerOccupiedSector;

    // Record data
    vector<vector<PData> > positionRecord;
    bool recPos;

    // Timing data
    high_resolution_clock::time_point start_time, end_time;

    // List of stat functions to use
    vector<StatFunc> statFunctions;
    vector<vector<pair<RealType,RealType> > > statFunctionData;
    vector<string> statFunctionName;

    // Performance tracking
    vector<pair<RealType,RealType> > performanceRecord;
    bool recPerf;
    high_resolution_clock::time_point last_record;

  };

}
#endif // __DATA_RECORD_HPP__
