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
#include "../objects/ScalarField.hpp"
#include "StatFunc.hpp"
#include "StatPlot.hpp"

namespace GFlow {

  // Forward declaration to Integrator
  class Integrator;
  // Forward declaration to VelocityVerletIntegrator
  class VelocityVerletIntegrator;

  /*
   * @class DataRecord
   * This class generates statistics from simulation data
   *
   */
  class DataRecord {
  public:
    // Default constructor
    DataRecord();

    // Initialize
    void initialize();

    // Start and end timing
    void startTiming();
    void endTiming();
   
    // Set lastRecord time for performance timing
    void markTime();

    // Set the setup time
    void setSetupTime(RealType t) { setupTime = t; }

    // Get the elapsed time
    double getElapsedTime() const;

    // Get the ratio of simulated time to real time
    RealType getRatio() const;

    // (Potentially) record data related to simulation data
    void record(SimData*, RealType); 
    void record(VelocityVerletIntegrator*, RealType);

    // Get sectorization data at the end
    void getSectorizationData(Sectorization*);

    // Write data
    void writeData(SimData* = nullptr) const;

    // Write run summary data
    void writeRunSummary(SimData* = nullptr, Integrator* = nullptr) const;

    /*** Mutators ***/

    // Set delay
    void setDelay(RealType d) { delay = d; }

    // Set command line options
    void setCommand(int, char**);

    // Set options
    void setRecPos(bool b)     { recPos = b; }
    void setRecOption(int i)   { recOption = i; }
    void setRecPerf(bool b)    { recPerf = b; }
    void setRecMvRatio(bool b) { recMvRatio = b; }
    void setRecDt(bool b)      { recDt = b; }
    void setRecDelay(bool b)   { recDelay = b; }
    void setRecBulk(bool b)    { recBulk = b; }
    void setRecBulkOutline(bool b) { recBulkOutline = b; }
    void setRecDisplacementField(bool b) { recDisplacementField = b; }
    void setRecPressField(bool b) { recPressField = b; }
    void setTrackDisplacement(bool b) { trackDisplacement = b; }

    // Whether to center the data
    void setCenter(bool c) { center = c; }

    // Set the write directory
    void setWriteDirectory(string w) { writeDirectory = w; }

    // Set the (expected) run time
    void setRunTime(RealType t) { runTime = t; }

    // Set the actual simulated time
    void setActualTime(RealType t) { actualTime = t; }

    // Add a stat function
    void addStatFunction(StatFunc, string);

    // Add a stat plo
    void addStatPlot(StatPlot, RPair, int, string);

    // Set whether to record the movement ratio
    void setRecMoveRatio(bool b) { recMvRatio = b; }

    /*** Accessors ***/

    // Get the simulation run time
    RealType getRunTime() const { return runTime; }

    // Get the number of stat functions
    int getNumberOfStatFunctions() const { return statFunctionData.size(); }

    // Get the data recorded by the stat function
    vector<pair<RealType, RealType> > getStatFunctionData(int) const;
    
    // Get the name of the stat function
    string getStatFunctionName(int) const;

    // Get the number of stat plots
    int getNumberOfStatPlots() const { return statPlotData.size(); }

    // --- Could have the other matching functions for stat plot here ---

    /*** Exception classes ***/
    struct BadStatFunction {
      BadStatFunction(int i) : value(i) {};
      int value;
    };

  private:
    // Private helper functions
    void writeParticleData(std::ofstream&, SimData*) const;
    void writeParticleChecks(std::ofstream&, SimData*) const;
    void recordByNumber(SimData*, vector<PData>&) const;
    void recordByVerletList(SimData*, vector<PData>&) const;
    void recordByVelocity(SimData*, vector<PData>&) const;
    void getBulkData(SimData*, const Bounds&, vector<RealType>&, ScalarField&, vector<pair<vec2,vec2> >&, RealType=0.015, RealType=0.015, RealType=0.01, RealType=1.) const;    
    inline void unite(int*, int, int) const;
    inline void createOutline(int*, int, int, RealType, RealType, Bounds, vector<pair<vec2,vec2> >&) const;
    inline int getHead(int*, int) const;
    inline void writeDisplacementData(SimData*) const;
    inline void writeDisplacementField(SimData*) const;
    inline void getPressureData(SimData*, const Bounds&, ScalarField&, RealType=0.1) const;

    // The command that was used
    vector<string> command;

    // The name of the directory we should write to 
    string writeDirectory;

    // How long between recording data
    RealType delay;

    // What time we last recorded data
    RealType lastRecord;

    // How long it took to parse
    RealType setupTime;

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
    int recOption; // Color particles according to what

    // Bulk data
    bool recBulk, recBulkOutline;
    ScalarField bulkField;
    vector<vector<pair<vec2,vec2> > > bulkOutline;
    vector<vector<RealType> > bulkVolumes;

    // Displacement field
    bool recDisplacementField;
    
    // Timing data
    high_resolution_clock::time_point start_time, end_time;

    // List of stat functions to use
    vector<StatFunc> statFunctions;
    vector<vector<RPair> > statFunctionData;
    vector<string> statFunctionName;

    // List of stat plots to use
    vector<StatPlot> statPlots;
    vector<vector<RPair> > statPlotData;
    vector<RPair> statPlotBounds;
    vector<string> statPlotName;

    // Displacement tracking
    bool trackDisplacement;
    vector<vec2> initialPositions;

    // Pressure field tracking
    bool recPressField;
    ScalarField pressField;

    // Performance tracking
    bool recPerf;
    high_resolution_clock::time_point last_record;
    int statRecPerf; // Which entry in statFunctionData refers to performance record

    bool recMvRatio;
    int statRecMvRatio;

    bool recDt;
    int statRecDt;

    bool recDelay;
    int statRecDelay;

    bool center; // Whether to center on the largest object
  };

}
#endif // __DATA_RECORD_HPP__
