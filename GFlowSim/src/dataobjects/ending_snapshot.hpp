#ifndef __ENDING_SNAPSHOT_HPP__GFLOW__
#define __ENDING_SNAPSHOT_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class EndingSnapshot : public DataObject {
  public:
    //! @brief Constructor
    EndingSnapshot(GFlow*);

    //! @brief Store the initial positions of the particles.
    virtual void pre_integrate() override;

    virtual void post_integrate() override;

    //! @brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    // Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

  private:
    //! @brief Helper function that stores all the neccesary data in a vector.
    inline void store_data(vector<RealType>&, vector<DataType>&);

    //! @brief Store a particular type of data.
    inline void get_data(vector<RealType>&, DataType, int);

    //! @brief Contains all the relevant data for a time step.
    //!
    //! Each time step contains a length ( [number of particles] * [dataWidth] ) array 
    //! for the positions of the particles in [DIMENSIONS] dimensions, and any other data we want to store.
    vector<RealType> final_data; 

    //! @brief Initial positions of particle.
    vector<RealType> initial_data;

    //! @brief The amount of data we collect per particle.
    int dataWidth;

    vector<DataType> data_types;
  };

}
#endif // __ENDING_SNAPSHOT_HPP__GFLOW__