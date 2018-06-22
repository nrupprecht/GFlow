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
  class DataMaster : protected Base {
  public:
    // Constructor
    DataMaster(GFlow *);

    // Destructor
    ~DataMaster();
    
    // Add a data object - we are subsequently in charge of the data object
    void addDataObject(DataObject*);

    // Call the corresponding routeens of the managed data objects - data
    // objects will collect data during one or more of these routines
    virtual void pre_step();
    virtual void pre_exchange();
    virtual void pre_forces();
    virtual void post_forces();
    virtual void post_step();

    // Do a coordinated write to a directory. Returns true if all writes were successful
    bool writeToDirectory(string);

    // GFlow is a friend class
    friend class GFlow;

  private:
    // The data objects we are responsible for
    vector<DataObject*> dataObjects;

  };

}
#endif // __DATA_MASTER_HPP__GFLOW__