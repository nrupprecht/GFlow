#ifndef __STORE_DATA_HPP__GFLOW__
#define __STORE_DATA_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {
  
  class StoreData {
  public:
    StoreData();

    //! \brief Compute the data places.
    void initialize(shared_ptr<class SimData>);

    void set_vector_data(const vector<string>&);
    void set_magnitude_data(const vector<string>&);
    void set_scalar_data(const vector<string>&);
    void set_integer_data(const vector<string>&);

    void set_data_boundary(const Bounds&);

    const vector<string>& get_vector_data() const;
    const vector<string>& get_scalar_data() const;
    const vector<string>& get_integer_data() const;

    //! \brief Store data to a vector.
    void store(vector<float>&);

    bool write(const string&, const vector<vector<float> >&);

    bool write(const string&, const vector<float>&);

    //! \brief Get the bounds.
    const Bounds& getBounds() const;
    //! \brief Get the data width
    int getDataWidth() const;
    int getDimensions() const;
    int getNTypes() const;

  private:

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
    bool write_processor_info;

    //! \brief The bounds data.
    Bounds bounds;
    //! \brief The data width.
    int dataWidth;
    //! \brief The number of dimensions.
    int sim_dimensions;
    //! \brief The number of types.
    int nTypes;

    //! \brief A pointer to the simdata class the data comes from.
    shared_ptr<class SimData> simData;

  };

}

#endif // __STORE_DATA_HPP__GFLOW__
