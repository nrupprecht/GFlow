#ifndef __DATA_MASTER_HPP__GFLOW__
#define __DATA_MASTER_HPP__GFLOW__

#include "dataobject.hpp"

namespace GFlowSimulation {

  /*
  *  @class DataMaster
  *
  *  Contains and coordinates data objects. GFlow has a single DataMaster, which in turn
  *  can have many data objects it uses to gather data. Data master is responsible for 
  *  saving data from its data objects to a unified output file system.
  *
  */
  class DataMaster : public Base {
  public:
    // Constructor
    DataMaster(GFlow *);

    // Destructor
    ~DataMaster();
    
    // Add a data object - we are subsequently in charge of the data object
    void addDataObject(DataObject*);

    // Start a timer
    void startTimer();

    // End the timer and add the new time to the record
    void endTimer();

    // Call the corresponding routeens of the managed data objects - data
    // objects will collect data during one or more of these routines
    virtual void pre_integrate();
    virtual void pre_step();
    virtual void pre_exchange();
    virtual void pre_forces();
    virtual void post_forces();
    virtual void post_step();
    virtual void post_integrate();

    // Do a coordinated write to a directory. Returns true if all writes were successful
    bool writeToDirectory(string);

    // Reset the time - use e.g. after relaxation step
    void resetTimer();

    // Set start recording time
    void setStartRecTime(RealType);

    // Set the fps of all the data objects
    void setFPS(RealType);

    // Set the fps of particular data objects
    void setFPS(int, RealType);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // --- Helper functions

    // Write a summary of the run to a text file
    inline bool writeSummary(string);

    // Write particle data to a stream
    inline void writeParticleData(std::ostream&);

    // Run time (real time)
    RealType run_time;
    // When the timer started
    high_resolution_clock::time_point start_time;
    // Whether the timer is running
    bool timing; 
    // Time to start taking data
    RealType startRecTime;

    // The data objects we are responsible for
    vector<DataObject*> dataObjects;

  };

}
#endif // __DATA_MASTER_HPP__GFLOW__