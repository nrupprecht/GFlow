#ifndef __POSITION_DATA_HPP__GFLOW__
#define __POSITION_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"
#include "../compute/store_data.hpp"

namespace GFlowSimulation {

  /*
  *  @class PositionData
  *  Records the position data of all the objects in the system.
  *  Do this in the post-step phase
  *
  */
  class PositionData : public DataObject {
  public:
    //! \brief Constructor.
    PositionData(GFlow*);

    //! \brief Store the initial positions of the particles.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    // Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

  private:

    //! \brief  The time steps of when the data was gathered
    vector<float> timeStamps;

    //! \brief Contains all the relevant data for a time step.
    //!
    //! Each time step contains a length ( [number of particles] * [dataWidth] ) array 
    //! for the positions of the particles in [DIMENSIONS] dimensions, and any other data we want to store.
    vector<vector<float> > positions; 

    //! \brief Initial positions of particle.
    vector<float> initial_data;

    // Data names and places
    vector<string> vector_data_entries;
    vector<string> scalar_data_entries;
    vector<string> integer_data_entries;

    StoreData storeData;
  };

}
#endif // __POSITION_DATA_HPP__GFLOW__