#include "sectorizationdata.hpp"
// Other files
#include "sectorization.hpp"

namespace GFlowSimulation {

  SectorizationData::SectorizationData(GFlow *gflow) : DataObject(gflow, "SecData") {};

  void SectorizationData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // Get the data
    auto &sectors = Base::sectorization->sectors;
    // Get the occupation numbers
    vector<int> occupation;
    for (int i=0; i<sectors.total(); ++i) occupation.push_back(sectors[i].size());
    // Add to the record
    sectorOccupation.push_back(occupation);
  }

  bool SectorizationData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);

    if(!PrintingUtility::writeVectorVectorToFile(sectorOccupation, dirName + dataName + ".csv")) return false;

    ofstream fout(dirName+"dims.csv");
    if (fout.fail()) return false;
    fout << PrintingUtility::toStrVec(Base::sectorization->dims) << endl;
    fout.close();

    return true;
  }

}