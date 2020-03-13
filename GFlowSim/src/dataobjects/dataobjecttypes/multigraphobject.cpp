#include "multigraphobject.hpp"

namespace GFlowSimulation {

  MultiGraphObject::MultiGraphObject(GFlow *gflow, const string& name, int nd) 
  : DataObject(gflow, name, DataObjectType::MULTIGRAPH), MultiGraphData(nd) {};

  MultiGraphObject::MultiGraphObject(GFlow *gflow, const string& name, const string& ax, const string& ay, int nd) 
    : DataObject(gflow, name, DataObjectType::MULTIGRAPH), MultiGraphData(ax, ay, nd) {};

  void MultiGraphObject::pre_integrate() {
    // Call parent class.
    DataObject::pre_integrate();
    // Reset data
    resetData();
  }

  bool MultiGraphObject::writeToFile(string fileName, bool make_directory) {
    // The name of the directory for this data
    string dirName = _correctDirName(fileName);
    // Create a directory for all the data
    if (make_directory) mkdir(dirName.c_str(), 0777);
    string name = dirName+dataName+"-"+toStr(object_counter)+".csv";
    return MultiGraphData::write_to_file(name);
  }

}