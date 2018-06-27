#include "verletlistdata.hpp"
// Other files
#include "force.hpp"

namespace GFlowSimulation {

  VerletListData::VerletListData(GFlow *gflow) : DataObject(gflow, "VLData"), nForces(0) {};

  VerletListData::~VerletListData() {};

  void VerletListData::pre_integrate() {
    nForces = forcesPtr->size();
  }

  void VerletListData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // Get the data
    for (auto &f : *Base::forcesPtr)
    verletLists.push_back(f->getVerletList());

  }

  bool VerletListData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);

    for (uint iter=0; iter<verletLists.size(); iter += nForces) {
      string fileHead = dirName + dataName + toStr(iter/nForces);
      for (int sub_iter=0; sub_iter<nForces; ++sub_iter) {
        string sub_fileName = (fileHead + "_" + toStr(sub_iter) + ".csv");
        PrintingUtility::writeVerletListToDirectory(verletLists.at(iter+sub_iter), sub_fileName);
      }
    }

    return true;
  }

}
