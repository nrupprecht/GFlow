#ifndef __LOAD_DATA_HPP__GFLOW__
#define __LOAD_DATA_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  class LoadData {
  public:
    LoadData();
    LoadData(string);

    bool load(string);

    const vector<vector<double> >& getData() const;
    const vector<double>& getData(int) const;
    int getEntries(int) const;
    const Bounds& getBounds() const;
    int getDataWidth() const;
    int getDimensions() const;
    int getNTypes() const;

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

    //! @brief The particle data.
    vector<vector<double> > data;
    //! @brief The bounds data.
    Bounds bounds;
    //! @brief The data width.
    int dataWidth;
    //! @brief The number of dimensions.
    int dimensions;
    //! @brief The number of types.
    int nTypes;
  };

}
#endif // __LOAD_DATA_HPP__GFLOW__