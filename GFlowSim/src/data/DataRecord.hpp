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
#include "StatFunc.hpp"

namespace GFlow {

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
    
    // Get the elapsed time
    double getElapsedTime();

    // (Potentially) record data
    void record(SimData*, RealType);

    // Write data
    void writeData(string, SimData* = nullptr);

    // Add a stat function
    void addStatFunction(StatFunc, string);

    /*** Accessors ***/

    // Get the number of stat functions
    int getNumberOfStatFunctions() { return statFunctionData.size(); }

    // Get the data recorded by the stat function
    vector<pair<RealType, RealType> > getStatFunctionData(int);
    
    // Get the name of the stat function
    string getStatFunctionName(int);

    /*** Exception classes ***/
    struct BadStatFunction {
      BadStatFunction(int i) : value(i) {};
      int value;
    };

  private:
    // How long between recording data
    RealType delay;

    // What time we last recorded data
    RealType lastRecord;

    // How many iterations we have recorded
    int recIter;

    // Record data
    vector<vector<PData> > positionRecord;

    // Timing data
    high_resolution_clock::time_point start_time, end_time;

    // List of stat functions to use
    vector<StatFunc> statFunctions;
    vector<vector<pair<RealType,RealType> > > statFunctionData;
    vector<string> statFunctionName;

  };

}
#endif // __DATA_RECORD_HPP__
