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
    
    // (Potentially) record data
    void record(SimData*, RealType);

    // Write data
    void writeData(string);
    
  private:
    // How long between recording data
    RealType delay;

    // What time we last recorded data
    RealType lastRecord;

    // How many iterations we have recorded
    int recIter;

    // Record data
    vector<vector<PData> > positionRecord;

  };

}
#endif // __DATA_RECORD_HPP__
