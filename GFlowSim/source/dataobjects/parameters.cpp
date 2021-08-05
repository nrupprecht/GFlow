#include <dataobjects/parameters.hpp>

using namespace GFlowSimulation;

Parameters::Parameters(GFlow *gflow)
    : DataObject(gflow, "Parameters") {};

void Parameters::addRecord(const string& name, RealType value) {
  data.insert(pair<string, RealType>(name, value));
}

bool Parameters::writeToFile(string fileName, bool useName) {
  // The name of the directory for this data
  string dirName = fileName;
  if (*fileName.rbegin() == '/') { // Make sure there is a /
    dirName += dataName + "/";
  }
  else {
    dirName += ("/" + dataName + "/");
  }

  // Write the data
  // Create a directory for all the data
  mkdir(dirName.c_str(), 0777);
  ofstream fout(dirName + dataName + ".csv");
  if (fout.fail()) {
    return false;
  }
  for (const auto& d : data) {
    fout << d.first << "," << d.second << endl;
  }
  fout.close();

  // Return success
  return true;
}

bool Parameters::find(const string &name, RealType &value) const {
  auto it = data.find(name);
  if (it == data.end()) {
    return false;
  }
  value = it->second;
  return true;
}