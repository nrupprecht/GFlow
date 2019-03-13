#include "load_data.hpp"

namespace GFlowSimulation {

  LoadData::LoadData() : bounds(Bounds(2)) {};

  LoadData::LoadData(string loadName) : bounds(Bounds(2)) {
    load(loadName);
  }

  bool LoadData::load(string loadName) {
    // Open file
    std::ifstream fin(loadName);
    if (fin.fail()) return false;
    // Get the data width, dimensions, number of times sampled, and the number of types
    int samples;
    dataWidth  = getNextNumber<int>(fin);
    dimensions = getNextNumber<int>(fin);
    samples    = getNextNumber<int>(fin);
    nTypes     = getNextNumber<int>(fin);
    
    // Get the dimensions - mins first, then maxes
    bounds = Bounds(dimensions);
    for (int i=0; i<dimensions; ++i) 
      bounds.min[i] = getNextNumber<RealType>(fin);
    for (int i=0; i<dimensions; ++i)
      bounds.max[i] = getNextNumber<RealType>(fin);

    // Get the types of data to expect
    string str;
    int nd = getNextNumber<int>(fin);
    for (int i=0; i<nd; ++i) {
      str = getNextString(fin);
      vector_data.push_back(str);
    }
    nd = getNextNumber<int>(fin);
    for (int i=0; i<nd; ++i) {
      str = getNextString(fin);
      scalar_data.push_back(str);
    }
    nd = getNextNumber<int>(fin);
    for (int i=0; i<nd; ++i) {
      str = getNextString(fin);
      integer_data.push_back(str);
    }

    // The first number in each line is the number of particles to expect
    for (int iter=0; iter<samples; ++iter) {
      // Get the length of data to expect
      int data_length = getNextNumber<int>(fin);
      // Get this iter's data
      vector<float> pdata;
      for (int i=0; i<data_length; ++i) {
        RealType datum = getNextNumber<float>(fin);
        pdata.push_back(datum);
      }
      // Store this iter's data vector
      data.push_back(pdata);
    }

    // Close the file stream
    fin.close();
    // Return success
    return true;
  }

  const vector<vector<float> >& LoadData::getData() const {
    return data;
  }

  int LoadData::getEntries(int i) const {
    return data.at(i).size()/dataWidth;
  }

  const vector<float>& LoadData::getData(int i) const {
    return data.at(i);
  } 

  const Bounds& LoadData::getBounds() const {
    return bounds;
  }

  int LoadData::getDataWidth() const {
    return dataWidth;
  }

  int LoadData::getDimensions() const {
    return dimensions;
  }

  int LoadData::getNTypes() const {
    return nTypes;
  }

  const vector<string>& LoadData::get_vector_data() {
    return vector_data;
  }

  const vector<string>& LoadData::get_scalar_data() {
    return scalar_data;
  }

  const vector<string>& LoadData::get_integer_data() {
    return integer_data;
  }

  inline string LoadData::getNextString(std::ifstream& fin) const {
    char c;
    fin.get(c);
    string str;
    while (!fin.eof() && c!=',' && c!='\n' && c!='\r') {
      str += c;
      // Get the next character
      fin.get(c);
    }
    // Convert to a type T and return
    return str;
  }

}