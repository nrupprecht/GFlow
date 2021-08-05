#ifndef __DATABASE_HPP__GFLOW__
#define __DATABASE_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  //! \brief An entry represents a column (or two, is use_percent==true) of a database. It is a purely abstract base class.
  struct Entry {
    Entry(const string& name) : entry_name(name) {};

    virtual ~Entry() {}

    //! \brief Resize the number of rows in the entry.
    virtual void resize(const int) = 0;

    //! \brief The number of rows in the entry.
    virtual int size() const = 0;

    //! \brief Return a string representation that contains the column name(s) for this entry.
    virtual string print_name() = 0;

    //! \brief Return a string representation of the data contained in the i-th row of this entry (column).
    virtual string print(const int) = 0;

    //! \brief Get the name of the entry.
    const string& getName() { 
      return entry_name; 
    }

  protected:
    //! \brief The name of the entry.
    string entry_name;

    //! \brief If true, then the print function will print the value, and what % of the total value the value constitutes.
    bool use_percent = false;
  };

  //! \brief A concrete entry, representing a column of a database. This inherits from Entry, so we 
  //! can store a collection of colums that have a different datatype polymorphically.
  template<typename T>
  struct GenericEntry : public Entry {
    GenericEntry(const string& name) : Entry(name) {};

    virtual void resize(const int s) override {
      data.resize(s);
    }

    virtual int size() const override {
      return data.size();
    }

    virtual string print_name() override {
      if (use_percent) return entry_name + "," + entry_name + " - %";
      else return entry_name;
    }

    virtual string print(const int i) override {
      if (0<=i && i<data.size()) {
        stringstream stream;
        stream << data[i];
        if (use_percent) {
          stream << "," << data[i]/total_value*100 << "%";
        }
        string str;
        stream >> str;
        return str;
      }
      else if (use_percent) return ",";
      else return "";
    }

    //! \brief Get a value from the column.
    T& get(int i) {
      return data.at(i);
    }

    void setData(const vector<T>& new_data) {
      data = new_data;
    }

    virtual void setUsePercent(bool use_per, T total = T(1)) {
      use_percent = use_per;
      total_value = total;
    }

    //! \brief Allows access to the type of this entry.
    typedef T type;

  private:
    //! \brief The underlying column data.
    vector<T> data;

    //! \brief Total value, for use when use_percent is true.
    T total_value = T(1);
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
      auto ptr = make_shared<GenericEntry<T> >(name);
      data_storage.push_back(ptr);
      ptr->setData(data);
    }

    //! \brief Print the database to a csv file, returns whether the write was successful or not.
    bool printToCSV(const string&);

    //! \brief Clear all the data. Does not reset n_rows.
    void clear();

    //! \brief Return a pointer to the most recently added entry.
    shared_ptr<Entry> last() const;

    //! \brief Set the row name.
    void setRowName(const string&);

  private:
    //! \brief The name of the first row. Default name is "Row."
    string row_name = "Row";

    //! \brief The vector of entries - each element is a column in the database.
    vector<shared_ptr<Entry> > data_storage;

    //! \brief The number of rows in the database.
    int n_rows = 0;
  };

}
#endif // __DATABASE_HPP__GFLOW__