#ifndef __ENDING_SNAPSHOT_HPP__GFLOW__
#define __ENDING_SNAPSHOT_HPP__GFLOW__

#include "../base/dataobject.hpp"
#include "../compute/store_data.hpp"

namespace GFlowSimulation {

  /*
  *  \brief The data object records data about the particles after integration.
  *
  *  This object works similarly to the position data object.
  *  \see PositionData
  */
  class EndingSnapshot : public DataObject {
  public:
    //! \brief Constructor
    EndingSnapshot(GFlow*);

    //! \brief Store the initial positions of the particles.
    virtual void pre_integrate() override;

    virtual void post_integrate() override;

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    // Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

  private:

    //! \brief Contains all the relevant data for a time step.
    //!
    //! Each time step contains a length ( [number of particles] * [dataWidth] ) array 
    //! for the positions of the particles in [DIMENSIONS] dimensions, and any other data we want to store.
    vector<RealType> final_data; 

    // Data names and places
    vector<string> vector_data_entries;
    vector<string> scalar_data_entries;
    vector<string> integer_data_entries;

    //! \brief Object for storing data.
    StoreData storeData;
  };

}
#endif // __ENDING_SNAPSHOT_HPP__GFLOW__