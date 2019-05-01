#ifndef __PARAMETERS_HPP__GFLOW__
#define __PARAMETERS_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  /**
  *  \brief A different type of data object, 
  *
  */
  class Parameters : public DataObject {
  public:
    //! \brief Default constructor.
    Parameters(GFlow*);

    //! \brief Add data to the parameter record.
    void addRecord(string, RealType);

    //! \brief Write the data to a file.
    bool writeToFile(string, bool);

  private:
    //! \brief Records of parameters.
    vector<pair<string, RealType> > data;
  };

}

#endif // __PARAMETERS_HPP__GFLOW__