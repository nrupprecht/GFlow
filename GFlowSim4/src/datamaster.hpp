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

    // GFlow is a friend class
    friend class GFlow;

    // Add a data object - we are subsequently in charge of the data object
    void addDataObject(DataObject*);

    // Have all data objects collect their data
    void collect();

    // Do a coordinated write to a directory
    void writeToDirectory(string);

  private:
    // The data objects we are responsible for
    vector<DataObject*> dataObjects;

  };

}
#endif // __DATA_MASTER_HPP__GFLOW__