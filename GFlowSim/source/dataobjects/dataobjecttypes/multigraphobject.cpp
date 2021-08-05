#include "dataobjects/dataobjecttypes/multigraphobject.hpp"

using namespace GFlowSimulation;

MultiGraphObject::MultiGraphObject(GFlow *gflow, const std::string &name, int nd)
    : DataObject(gflow, name, DataObjectType::MULTIGRAPH), MultiGraphData(nd) {};

MultiGraphObject::MultiGraphObject(GFlow *gflow,
                                   const std::string &name,
                                   const std::string &ax,
                                   const std::string &ay, int nd)
    : DataObject(gflow, name, DataObjectType::MULTIGRAPH), MultiGraphData(ax, ay, nd) {};

void MultiGraphObject::pre_integrate() {
  // Call parent class.
  DataObject::pre_integrate();
  // Reset data
  resetData();
}

bool MultiGraphObject::writeToFile(std::string fileName, bool useName) {
  // The name of the directory for this data
  string dirName = _correctDirName(fileName);
  mkdir(dirName.c_str(), 0777);
  return MultiGraphData::write_to_file(dirName, dataName + "-" + toStr(object_counter));
}
