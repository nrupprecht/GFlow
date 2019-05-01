#include "parameters.hpp"

namespace GFlowSimulation {

  Parameters::Parameters(GFlow *gflow) : DataObject(gflow, "Parameters") {};

  void Parameters::addRecord(string name, RealType value) {
    data.push_back(pair<string, RealType>(name, value));
  }

  bool Parameters::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    // Write the data
    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);
    ofstream fout(dirName+dataName+".csv");
    if (fout.fail()) return false;
    for (auto d : data) {
      fout << d.first << "," << d.second << endl;
    }
    fout.close();

    // Return success
    return true;
  }

}