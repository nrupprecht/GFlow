#ifndef __PARTICLE_STORE_DATA_HPP__GFLOW__
#define __PARTICLE_STORE_DATA_HPP__GFLOW__

#include "../../../gflow.hpp"
#include "../../../compute/store_data.hpp"

namespace GFlowSimulation {

  /** \brief A data base for objects that need to collect particle data.
   *  
   *  Data base objects are used so that if I need a data object to inherit from multiple "data object"
   *  like functionalities, I don't get triangle inheritance.
   * 
   */
  class ParticleStoreData {
  public:
    //! \brief Request to store a vector data entry.
    void add_vector_data_entry(string);

    //! \brief Request to store the magnitude of a vector data entry.
    void add_magnitude_data_entry(string);

    //! \brief Request to store a scalar data entry.
    void add_scalar_data_entry(string);

    //! \brief Request to store an integer data entry.
    void add_integer_data_entry(string);

    //! \brief Clear all the requested data entries, of every type.
    void clear_all_data_entries();

  protected:

    //! \brief Vector data names.
    vector<string> vector_data_entries;

    //! \brief Vector magnitude data names.
    vector<string> magnitude_data_entries;

    //! \brief Scalar data names.
    vector<string> scalar_data_entries;

    //! \brief Integer data names.
    vector<string> integer_data_entries;

    //! \brief A store data object.
    StoreData storeData;

    //! \brief A function that can be used to only record some particles.
    std::function<bool(shared_ptr<SimData>, int)> select_function;

  };

}

#endif // __PARTICLE_STORE_DATA_HPP__GFLOW__