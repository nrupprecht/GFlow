#include <dataobjects/averagedata.hpp>

using namespace GFlowSimulation;

AverageData::AverageData(GFlow *gflow)
    : DataObject(gflow, "Averages"), dataWidth(4) {};

AverageData::~AverageData() {
  for (auto &v : data) {
    delete v;
  }
  data.clear();
}

void AverageData::post_step() {
  // Only record if enough time has gone by
  if (!DataObject::_check()) {
    return;
  }
  // Store data
  RealType time = gflow->getElapsedTime();
  auto dat = new RealType[dataWidth];
  // Put in values
  dat[0] = time;
  dat[1] = 0;
  dat[2] = 0;
  dat[3] = 0;
  // Get and store data
  int size = simData->size_owned(), number = simData->number_owned();
  for (int n = 0; n < size; ++n) {
    if (simData->Type(n) < 0) {
      continue;
    }
    dat[1] += magnitudeVec(simData->X(n), sim_dimensions);
    dat[2] += magnitudeVec(simData->V(n), sim_dimensions);
    dat[3] += magnitudeVec(simData->F(n), sim_dimensions);
  }
  // Put in values
  dat[1] /= number;
  dat[2] /= number;
  dat[3] /= number;
  // Store
  data.push_back(dat);
}

bool AverageData::writeToFile(string fileName, bool useName) {
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
  ofstream fout(dirName + dataName + ".csv");
  if (fout.fail()) {
    return false;
  }
  for (auto d : data) {
    for (int i = 0; i < dataWidth; ++i) {
      fout << d[i];
      if (i != dataWidth - 1) {
        fout << ",";
      }
    }
    fout << endl;
  }
  fout.close();

  // Return success
  return true;
}
