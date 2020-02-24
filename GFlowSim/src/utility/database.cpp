#include "database.hpp"

namespace GFlowSimulation {

  DataBase::DataBase(const int nr) : n_rows(nr) {};

  //! \brief Set the number of rows in the database.
  void DataBase::setNRows(const int nr) {

  }

  //! \brief Add a column to the database.
  void DataBase::addColumn(const string& name, const char type) {
    if (type=='I') {
      auto ptr = make_shared<GenericEntry<int> >(name);
      ptr->resize(n_rows);
      data_storage.push_back(ptr);
    }
    else if (type=='R') {
      auto ptr = make_shared<GenericEntry<real> >(name);
      ptr->resize(n_rows);
      data_storage.push_back(ptr);
    }
    else cout << "Error - unrecognized type character.\n";
  }

  //! \brief Print the database to a csv file.
  void DataBase::printToCSV(const string& fileName) {

  }

}