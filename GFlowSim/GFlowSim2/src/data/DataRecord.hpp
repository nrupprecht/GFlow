/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __DATA_RECORD_HPP__
#define __DATA_RECORD_HPP__

// Includes

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
    
    // Update the data record's time
    void update(RealType dt) { delayTimer += dt; }

    // (Potentially) record data
    void record(SimData*);
    
  private:
    // How long between recording data
    RealType delay;
    // How long since data was last recorded
    RealType delayTimer;

  };

}
#endif // __DATA_RECORD_HPP__
