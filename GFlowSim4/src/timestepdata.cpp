#include "timestepdata.hpp"

namespace GFlowSimulation {

  TimeStepData::TimeStepData(GFlow *gflow) : DataObject(gflow, "Timestep") {};

   void TimeStepData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Store data
    RealType time = Base::gflow->getElapsedTime();
    data.push_back(RPair(time, gflow->getDT()));
  }

  bool TimeStepData::writeToFile(string fileName, bool useName) {
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
    for (auto ts : data)
      fout << ts.first << "," << ts.second << endl;
    fout.close();

    // Return success
    return true;
  }

}