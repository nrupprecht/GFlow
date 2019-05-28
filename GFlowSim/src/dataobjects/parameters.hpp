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

    //! \brief Try to get a parameter by name. Returns true if successful.
    bool find(const string&, RealType&) const;

  private:
    //! \brief Records of parameters.
    std::map<string, RealType> data;
  };

}

#endif // __PARAMETERS_HPP__GFLOW__