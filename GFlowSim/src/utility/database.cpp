#include "database.hpp"

namespace GFlowSimulation {

  DataBase::DataBase(const int nr) : n_rows(nr) {};

  //! \brief Set the number of rows in the database.
  void DataBase::setNRows(const int nr) {
    n_rows = nr;
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
  bool DataBase::printToCSV(const string& fileName) {

    ofstream fout(fileName);
    if (fout.fail()) {
      std::cerr << "Failed to open file in DataBase::printToCSV.\n";
      return false;
    }
  
    // Print labels, then print each row.

    fout << row_name;
    for (auto ptr : data_storage) { 
      fout << "," << ptr->print_name();
    }
    fout << "\n";

    for (int i=0; i<n_rows; ++i) {    
      fout << i;
      for (auto ptr : data_storage) {
        fout << "," << ptr->print(i);
      }
      fout << "\n";
    }

    fout.close();
    return true;
  }

  void DataBase::clear() {
    data_storage.clear();
  }

  shared_ptr<Entry> DataBase::last() const {
    if (!data_storage.empty()) {
      return data_storage[data_storage.size()-1];
    }
    else return nullptr;
  }

  void DataBase::setRowName(const string& name) {
    row_name = name;
  }

}