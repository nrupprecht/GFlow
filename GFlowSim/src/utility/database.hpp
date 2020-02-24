#ifndef __DATABASE_HPP__GFLOW__
#define __DATABASE_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  //! \brief An entry represents a column of a database. It is a purely abstract class.
  struct Entry {
    Entry(const string& name) : entry_name(name) {};

    //! \brief Resize the number of rows in the entry.
    virtual void resize(const int s)=0;

    //! \brief Get the name of the entry.
    const string& getName() { 
      return entry_name; 
    }

  private:
    //! \brief The name of the entry.
    string entry_name;
  };

  //! \brief A concrete entry, representing a column of a database. This inherits from Entry, so we 
  //! can store a collection of colums that have a different datatype polymorphically.
  template<typename T>
  struct GenericEntry : public Entry {
    GenericEntry(const string& name) : Entry(name) {};

    virtual void resize(const int s) override {
      data.resize(s);
    }

    //! \brief Get a value from the column.
    T& get(int i) {
      return data.at(i);
    }

    //! \brief Allows access to the type of this entry.
    typedef T type;

  private:
    //! \brief The underlying column data.
    vector<T> data;
  };



  class DataBase {
  public:
    DataBase() = default;

    DataBase(const int nr);

    //! \brief Set the number of rows in the database.
    void setNRows(const int);

    //! \brief Add a column to the database.
    void addColumn(const string&, const char);

    template<typename T>
    void addColumn(const string& name, const vector<T>& data) {
      data_storage.push_back(make_shared<GenericEntry<T> >(name, n_rows));
    }

    //! \brief Print the database to a csv file.
    void printToCSV(const string&);

  private:
    //! \brief The vector of entries - each element is a column in the database.
    vector<shared_ptr<Entry> > data_storage;

    //! \brief The number of rows in the database.
    int n_rows = 0;
  };

}
#endif // __DATABASE_HPP__GFLOW__