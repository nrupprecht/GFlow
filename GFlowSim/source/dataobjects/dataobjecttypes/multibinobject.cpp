#include "dataobjects/dataobjecttypes/multibinobject.hpp"
// Other files.
#include <utility>

using namespace GFlowSimulation;

MultiBinObject::MultiBinObject(GFlow *gflow, string name, int nd)
    : DataObject(gflow, std::move(name), DataObjectType::MULTIBIN), ndata(max(nd, 1)) {
  // Write data flags are true by default.
  write_data = vector<bool>(nd, true);
  // Default y axis name.
  axis_y = vector<string>(nd, "y");
  // Set up data array.
  resetData();
}

void MultiBinObject::pre_integrate() {
  resetData();
}

bool MultiBinObject::writeToFile(string fileName, bool useName) {
  // Check if there's anything to do
  if (ndata_points == 0) {
    return true;
  }
  // The name of the directory for this data
  string dirName = _correctDirName(fileName);
  // Create a directory for all the data
  mkdir(dirName.c_str(), 0777);
  ofstream fout(dirName + dataName + "-" + toStr(object_counter) + ".csv");
  if (fout.fail()) {
    return false;
  }
  // Print out the data
  for (int i = 0; i < counts.size(); ++i) {
    fout << labels[i];
    for (int j = 0; j < ndata; ++j) {
      if (j == 0 || write_data[j - 1]) {
        // Prepend with a "," (if not the first entry) so there will not be null entries if the last entry should not be printed.
        fout << ",";
        fout << counts[i][j];
      }
    }
    fout << endl;
  }
  fout.close();
  // Print out all the axes labels.
  fout.open(dirName + "axes.csv");
  if (fout.fail()) {
    return false;
  }
  for (int i = 0; i < ndata; ++i) {
    if (write_data[i]) {
      fout << axis_x << "," << axis_y[i] << endl;
    }
  }
  fout.close();

  // Return success
  return true;
}

void MultiBinObject::setNBins(int nb) {
  if (nb != nbins) {
    nbins = nb;
    resetData();
    labelBins(min_label, max_label);
  }
}

void MultiBinObject::increment(int bin, int entry, bool inc) {
  ++counts[bin][entry];
  if (inc) {
    ++ndata_points;
  }
}

void MultiBinObject::increment(int bin, bool inc) {
  ++counts[bin][0];
  if (inc) {
    ++ndata_points;
  }
}

void MultiBinObject::resetData() {
  counts = vector<vector<long> >(nbins, vector<long>(ndata));
  labels = vector<RealType>(nbins);
  ndata_points = 0;
}

void MultiBinObject::labelBins(RealType min, RealType max) {
  // If the vector is not the right size, resize it.
  if (labels.size() != nbins) {
    labels = vector<RealType>(nbins);
  }
  // Create labels.
  RealType dx = (max - min) / (nbins - 1);
  for (int i = 0; i < nbins; ++i) {
    labels[i] = min + i * dx;
  }
  // Set min/max.
  min_label = min;
  max_label = max;
}
