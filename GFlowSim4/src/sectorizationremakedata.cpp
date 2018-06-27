#include "sectorizationremakedata.hpp"
// Other files
#include "sectorization.hpp"

namespace GFlowSimulation {
  
  SectorizationRemakeData::SectorizationRemakeData(GFlow *gflow) : DataObject(gflow, "SecRemake") {};

  void SectorizationRemakeData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    int remakes = Base::sectorization->getNumberOfRemakes();
    RealType time = Base::gflow->getElapsedTime();
    remakeData.push_back(RIPair(time, remakes));
  }

  bool SectorizationRemakeData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = _correctDirName(fileName);

    // Create a directory for all the data
    _makeDir(dirName);

    // Write data
    ofstream fout(dirName+dataName+".csv");
    if (fout.fail()) return false;
    for (auto c : remakeData)
      fout << c.first << "," << c.second << endl;
    fout.close();
    // Return success
    return true;
  }

}