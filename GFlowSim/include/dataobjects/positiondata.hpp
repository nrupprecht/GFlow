#ifndef __POSITION_DATA_HPP__GFLOW__
#define __POSITION_DATA_HPP__GFLOW__

#include "../base/dataobject.hpp"
//#include "../compute/store_data.hpp"

#include "dataobjecttypes/data-base/particle-store-data.hpp"

namespace GFlowSimulation {

  /*
  *  \class PositionData
  *  Records the position data of all the objects in the system.
  *  Do this in the post-step phase.
  *
  */
  class PositionData : public DataObject, public ParticleStoreData {
  public:
    //! \brief Constructor.
    PositionData(GFlow*);

    //! \brief Store the initial positions of the particles.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    //!
    //! Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

    // //! \brief Request to store a vector data entry.
    // void add_vector_data_entry(string);
    // //! \brief Request to store the magnitude of a vector data entry.
    // void add_magnitude_data_entry(string);
    // //! \brief Request to store a scalar data entry.
    // void add_scalar_data_entry(string);
    // //! \brief Request to store an integer data entry.
    // void add_integer_data_entry(string);

    // //! \brief Clear all the requested data entries, of every type.
    // void clear_all_data_entries();

    // The demon modifer is a friend. This way, the demon can correct the video 
    // when it resets the clock and closes the door.
    friend class Demon;

  protected:

    //! \brief The time steps of when the data was gathered
    vector<float> timeStamps;

    //! \brief Contains all the relevant data for a time step.
    //!
    //! Each time step contains a length ( [number of particles] * [dataWidth] ) array 
    //! for the positions of the particles in [DIMENSIONS] dimensions, and any other data we want to store.
    vector<vector<float> > positions; 

    //! \brief Initial positions of particle.
    vector<float> initial_data;

    // //! \brief Vector data names.
    // vector<string> vector_data_entries;
    // //! \brief Vector magnitude data names.
    // vector<string> magnitude_data_entries;
    // //! \brief Scalar data names.
    // vector<string> scalar_data_entries;
    // //! \brief Integer data names.
    // vector<string> integer_data_entries;

    // //! \brief A store data object.
    // StoreData storeData;

    // //! \brief A function that can be used to only record some particles.
    // std::function<bool(shared_ptr<SimData>, int)> select_function;
  };

}
#endif // __POSITION_DATA_HPP__GFLOW__
