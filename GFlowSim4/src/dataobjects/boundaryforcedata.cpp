#include "boundaryforcedata.hpp"

namespace GFlowSimulation {

  BoundaryForceData::BoundaryForceData(GFlow *gflow) : DataObject(gflow, "BDForce") {};

  void BoundaryForceData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Store data
    bForces.push_back(RPair( Base::gflow->getElapsedTime(), 
                            Base::gflow->getBoundaryForce()));
  }

  bool BoundaryForceData::writeToFile(string fileName, bool useName) {
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
    for (const auto bf : bForces)
      fout << bf.first << "," << bf.second << endl;
    fout.close();

    // Return success
    return true;
  }

  RealType BoundaryForceData::getAverage() const {
    RealType force = 0;
    for (const auto bf : bForces)
      force += bf.second;
    return force/bForces.size();
  }

}