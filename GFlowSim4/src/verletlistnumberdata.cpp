#include "verletlistnumberdata.hpp"
// Other files
#include "force.hpp"

namespace GFlowSimulation {

  VerletListNumberData::VerletListNumberData(GFlow *gflow) : DataObject(gflow, "VLNumbers") {};

  void VerletListNumberData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    // Get the data
    RealType time = Base::gflow->getElapsedTime();
    int nheads(0), nverlet(0);
    for (auto &f : *Base::forcesPtr) {
      nheads += f->vlHSize();
      nverlet += f->vlSize();
    }
    verletNumbers.push_back(pair<RealType, IPair>(time, IPair(nheads, nverlet)));
  }

  bool VerletListNumberData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);

    ofstream fout(dirName+dataName+".csv");
    if (fout.fail()) return false;
    for (auto d : verletNumbers)
      fout << d.first << "," << d.second.first << "," << d.second.second << endl;
    fout.close();
    // Return success
    return true;
  }

}