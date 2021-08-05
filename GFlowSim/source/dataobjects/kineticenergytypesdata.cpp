#include <dataobjects/kineticenergytypesdata.hpp>
// Other files
#include <base/simdata.hpp>
#include <utility/vectormath.hpp>

using namespace GFlowSimulation;

KineticEnergyTypesData::KineticEnergyTypesData(GFlow *gflow, bool ave)
    : DataObject(gflow, "KETypes"), useAve(ave) {
  ntypes = gflow->getNTypes();
  keData = new vector<RPair>[ntypes];
};

KineticEnergyTypesData::KineticEnergyTypesData(GFlow *gflow, int types, bool ave)
    : DataObject(gflow, "KETypes"), ntypes(types), useAve(ave) {
  keData = new vector<RPair>[ntypes];
};

KineticEnergyTypesData::~KineticEnergyTypesData() {
  if (keData) {
    delete[] keData;
  }
}

void KineticEnergyTypesData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }
  RealType time = gflow->getElapsedTime();
  // Create an entry
  for (int ty = 0; ty < ntypes; ++ty) {
    keData[ty].push_back(RPair(time, 0));
  }
  // Get and store data
  auto v = simData->V();
  auto im = simData->Im();
  int size = simData->size_owned();
  for (int n = 0; n < size; ++n) {
    int type = simData->Type(n);
    if (-1 < type && type < ntypes && im[n] > 0) {
      keData[type].rbegin()->second += sqr(v(n), sim_dimensions) / im[n];
    }
  }
  // Multiply by 1/2
  for (int ty = 0; ty < ntypes; ++ty) {
    keData[ty].rbegin()->second *= 0.5;
    if (useAve) {
      keData[ty].rbegin()->second /= simData->number();
    }
  }
}

bool KineticEnergyTypesData::writeToFile(string fileName, bool useName) {
  // The name of the directory for this data
  string dirName = fileName;
  if (*fileName.rbegin() == '/') { // Make sure there is a /
    dirName += dataName + "/";
  }
  else {
    dirName += ("/" + dataName + "/");
  }

  // Write the data
  // Create a directory for all the data
  mkdir(dirName.c_str(), 0777);
  for (int ty = 0; ty < ntypes; ++ty) {
    ofstream fout(dirName + dataName + toStr(ty) + ".csv");
    if (fout.fail()) {
      return false;
    }
    for (auto ke : keData[ty]) {
      fout << ke.first << "," << ke.second << endl;
    }
    fout.close();
  }

  // Return success
  return true;
}
