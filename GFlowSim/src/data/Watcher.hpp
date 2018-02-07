/*
 * Author: Nathaniel Rupprecht
 * Start Data: Feb, 6, 2018
 *
 */

#ifndef __WATCHER_HPP__
#define __WATCHER_HPP__

#include "../../include/Utility.hpp"

namespace GFlow {

  // Forward declarations
  class DataRecord;
  class SimData;

  class Watcher {
  public:
    // Constructor
    Watcher(int, bool, DataRecord*, string);
    // Set data record
    void setDataRecord(DataRecord*);
    // Record data 
    virtual void record() = 0;
  protected:
    // Get the sim data from dataRecord
    SimData* getSimData();
    // Particle id
    int id;
    // True for particle, false for wall
    bool type;
    // Pointer to the data record
    DataRecord *dataRecord;
    // The slot in dataRecord to use
    int slot;
    // The name of this stat
    string name;
  };

  class NetForceWatcher : public Watcher {
  public:
    // Constructor
    NetForceWatcher(int, bool, DataRecord* = nullptr);
    // Record Data
    virtual void record();
  };

  class LKEWatcher : public Watcher {
  public:
    // Constructor
    LKEWatcher(int, bool, DataRecord* = nullptr);
    // Record Data
    virtual void record();
  };

  class PyWatcher : public Watcher {
  public:
    // Constructor
    PyWatcher(int, bool, DataRecord* = nullptr);
    // Record Data
    virtual void record();
  };

}

#endif // __WATCHER_HPP__
