#ifndef __LOAD_DATA_HPP__GFLOW__
#define __LOAD_DATA_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  class LoadData {
  public:
    LoadData();
    LoadData(const string&);

    bool load(const string&);

    const vector<vector<float> >& getData() const;
    const vector<float>& getData(int) const;
    int getEntries(int) const;
    const Bounds& getBounds() const;
    int getDataWidth() const;
    int getDimensions() const;
    int getNTypes() const;

    const vector<string>& get_vector_data();
    const vector<string>& get_scalar_data();
    const vector<string>& get_integer_data();

  private:

    template<typename T> inline T getNextNumber(std::ifstream& fin) const {
      char c;
      fin.get(c);
      stringstream stream;
      while (!fin.eof() && c!=',' && c!='\n') {
        if (isdigit(c) || c=='.' || c=='-' || c=='e' || c=='E') stream << c;
        else;
        // Get the next character
        fin.get(c);
      }
      // Convert to a type T and return
      T val;
      stream >> val;
      return val;
    }

    inline string getNextString(std::ifstream& fin) const;

    //! @brief The particle data.
    vector<vector<float> > data;
    //! @brief The bounds data.
    Bounds bounds;
    //! @brief The data width.
    int dataWidth;
    //! @brief The number of dimensions.
    int dimensions;
    //! @brief The number of types.
    int nTypes;

    // The names of the data that we are reading in.
    vector<string> vector_data;
    vector<string> scalar_data;
    vector<string> integer_data;
  };

}
#endif // __LOAD_DATA_HPP__GFLOW__