#ifndef __STORE_DATA_HPP__GFLOW__
#define __STORE_DATA_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {
  
  class StoreData {
  public:

    //! \brief Compute the data places.
    void initialize(const shared_ptr<class SimData>&);

    //! \brief Set the list of vector data elements to store.
    void set_vector_data(const vector<string>&);
    //! \brief Set the list of vector magnitudes to store.
    void set_magnitude_data(const vector<string>&);
    //! \brief Set the list of scalar data elements to store.
    void set_scalar_data(const vector<string>&);
    //! \brief Set the list of integer data elements to store.
    void set_integer_data(const vector<string>&);

    //! \brief Set the bounds for which particle data should be captured.
    void set_data_boundary(const Bounds&);

    const vector<string>& get_vector_data() const;
    const vector<string>& get_scalar_data() const;
    const vector<string>& get_integer_data() const;

    //! \brief Store data to a vector.
    void store(vector<float>&);
    //! \brief Store data to a vector of particles that fit the criteria of some boolean function.
    void store(vector<float>&, const std::function<bool(std::shared_ptr<SimData>, int)>&);

    //! \brief Write a collection of data frames to a csv file.
    bool write(const string&, const vector<vector<float> >&);
    //! \brief Write a single data frame to a csv file.
    bool write(const string&, const vector<float>&);

    //! \brief Get the bounds.
    const Bounds& getBounds() const;
    //! \brief Get the data width
    int getDataWidth() const;
    //! \brief Get the dimensionality.
    int getDimensions() const;
    //! \brief Get the number of types.
    int getNTypes() const;

  private:
    //! \brief Helper function that writes information about the types of data that the csv file contains.
    inline void writeHeader(ofstream&, int) const;

    // Data names and places
    vector<string> vector_data_entries;
    vector<string> magnitude_data_entries;
    vector<string> scalar_data_entries;
    vector<string> integer_data_entries;

    vector<int> vector_data_positions;
    vector<int> magnitude_data_positions;
    vector<int> scalar_data_positions;
    vector<int> integer_data_positions;

    //! \brief Whether to write processor ownership related info. This would be the last entry.
    bool write_processor_info = true;

    //! \brief The bounds data.
    Bounds bounds = Bounds(2);
    //! \brief The data width.
    int dataWidth = 0;
    //! \brief The number of dimensions.
    int sim_dimensions = 0;
    //! \brief The number of types.
    int nTypes = 0;

    //! \brief A pointer to the simdata class the data comes from.
    shared_ptr<class SimData> simData;

  };

}

#endif // __STORE_DATA_HPP__GFLOW__
