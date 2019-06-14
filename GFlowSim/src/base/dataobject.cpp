#include "dataobject.hpp"

namespace GFlowSimulation {

  DataObject::DataObject(GFlow *gflow, const string& name) 
    : Base(gflow), dataName(name), delay(1./20.), lastRecording(-10.), type(DataObjectType::GENERAL), locals_changed(false) {
      object_counter = 1;
    };

  DataObject::DataObject(GFlow *gflow, const string& name, DataObjectType t)
    : Base(gflow), dataName(name), delay(1./20.), lastRecording(-10.), type(t), locals_changed(false) {
      object_counter = 1;
    };

  void DataObject::pre_integrate() {
    lastRecording = -10;
  }

  const string& DataObject::getName() const {
    return dataName;
  }

  DataObjectType DataObject::getType() const {
    return type;
  }

  RealType DataObject::getLastRecording() const {
    return lastRecording;
  }

  void DataObject::setFPS(RealType fps) {
    delay = 1./fps;
  }

  void DataObject::setLocalsChanged(bool r) {
    locals_changed = r;
  }

  int DataObject::getObjectCounter() {
    return object_counter;
  }

  void DataObject::setDataName(const string& d) {
    dataName = d;
  }

  void DataObject::setLastRecording(RealType t) {
    lastRecording = t;
  }

  string DataObject::_correctDirName(string fileName) {
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"-"+toStr(object_counter)+"/");
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
