#include "dataobject.hpp"

namespace GFlowSimulation {

  DataObject::DataObject(GFlow *gflow, string name) : Base(gflow), dataName(name), delay(1./20.), lastRecording(-10.) {};

  string DataObject::getName() {
    return dataName;
  }

  void DataObject::setFPS(RealType fps) {
    delay = 1./fps;
  }

  string DataObject::_correctDirName(string fileName) {
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");
    return dirName;
  }

  void DataObject::_makeDir(string dirName) {
    mkdir(dirName.c_str(), 0777);
  }

  bool DataObject::_check() {
    // Only record if enough time has gone by
    RealType time = Base::gflow->getElapsedTime();
    // If not enough time has gone by, return false
    if (time-lastRecording<delay) return false;
    // Otherwise, set last recording, return true
    lastRecording = time;
    return true;
  }

}
