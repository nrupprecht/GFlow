#include "dataobjects/dataobjecttypes/volumeplotobject2d.hpp"
// Other files.
#include "utility/vectormath.hpp"

using namespace GFlowSimulation;

// --- VolumePlotObject2D ---

VolumePlotObject2D::VolumePlotObject2D(GFlow *gflow, const string &name, int width)
    : DataObject(gflow, name, DataObjectType::VOLUMEPLOT),
      data_width(width) {
  if (width <= 0) {
    throw Exception("Data width is <= 0 in VolumePlotObject2D constructor.");
  }
  // Create some default entry names.
  for (int i = 0; i < data_width; ++i) {
    entry_names.push_back("value-" + toStr(i));
  }
};

void VolumePlotObject2D::pre_integrate() {
  // If needed, set the focus bounds to be the simulation bounds.
  if (gather_bounds.vol() == 0) {
    gather_bounds = gflow->getBounds();
  }

  // If binX, binY not set, choose automatically.
  int new_bin_x = binX, new_bin_y = binY;
  if (binX == 0 || binY == 0) {
    RealType wx = gather_bounds.wd(0);
    RealType wy = gather_bounds.wd(1);
    // Set up bins so the larger direction has default_max_bins # of bins.
    if (wx < wy) {
      new_bin_y = default_max_bins;
      RealType bin_width = wy / new_bin_y;
      new_bin_x = static_cast<int>(wx / bin_width);
    }
    else { // wy <= wx
      new_bin_x = default_max_bins;
      RealType bin_width = wx / new_bin_x;
      new_bin_y = static_cast<int>(wy / bin_width);
    }
  }

  // Set bins
  setBins(new_bin_x, new_bin_y);

  // Reset.
  recorded_frames = 0;
}

bool VolumePlotObject2D::writeToFile(string fileName, bool useName) {
  // Check if there's anything to do
  if (binning.empty() || counts.empty() || data_width <= 0) {
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
  RealType dx = gather_bounds.wd(0) / binX, dy = gather_bounds.wd(1) / binY;
  for (int bx = 0; bx < binX; ++bx) {
    for (int by = 0; by < binY; ++by) {
      // Print coordinates of the center of the bin.
      fout << gather_bounds.min[0] + (bx + 0.5) * dx << "," << gather_bounds.min[1] + (by + 0.5) * dy;
      // Print out data
      for (int d = 0; d < data_width; ++d) {
        fout << "," << (counts[bx][by] > 0 ? binning[bx][by][d] / counts[bx][by] : 0);
      }
      // Finally, print out counts, if requested.
      if (print_counts) {
        fout << "," << (recorded_frames > 0 ? static_cast<float>(counts[bx][by]) / recorded_frames : counts[bx][by])
             << endl;
      }
    }
  }
  fout.close();

  // Open file for printing the entry names.
  fout.open(dirName + "entrynames.csv");
  if (fout.fail()) {
    return false;
  }
  for (int i = 0; i < data_width; ++i) {
    fout << entry_names[i];
    if (i != data_width - 1) {
      fout << ",";
    }
  }
  fout.close();

  // Return success
  return true;
}

void VolumePlotObject2D::setBins(int bx, int by) {
  // We need to resize.
  if (binX != bx || binY != by) {
    binning = vector<vector<Vec> >(bx, vector<Vec>(by, Vec(data_width)));
    counts = vector<vector<int> >(bx, vector<int>(by, 0));
  }
    // Don't need to resize. Just set all values to zero.
  else {
    // Set binning to zero.
    for (auto &row : binning) {
      for (auto &entry : row) {
        entry.zero();
      }
    }
    // Set counts to zero.
    for (auto &row : counts) {
      for (auto &entry : row) {
        entry = 0;
      }
    }
  }
  // Set binX, binY
  binX = by;
  binY = by;
  // Calculate invdx, invdy
  invdx = binX / gather_bounds.wd(0);
  invdy = binY / gather_bounds.wd(1);
}

void VolumePlotObject2D::setPrintCounts(bool p) {
  print_counts = p;
}

void VolumePlotObject2D::addToBin(int x, int y, int d, RealType v, bool inc) {
  binning[x][y][d] += v;
  if (inc) {
    ++counts[x][y];
  }
}

void VolumePlotObject2D::addToBin(RealType px, RealType py, int d, RealType v, bool inc) {
  // Calculate bin index.
  int x = static_cast<int>(invdx * (px - gather_bounds.min[0]));
  int y = static_cast<int>(invdy * (py - gather_bounds.min[1]));
  // Add to bin by bin index.
  if (0 <= x && x < binX && 0 <= y && y < binY) {
    addToBin(x, y, d, v, inc);
  }
}

void VolumePlotObject2D::addToBin(RealType px, RealType py, Vec &v) {
  // Calculate bin index.
  int x = static_cast<int>(invdx * (px - gather_bounds.min[0]));
  int y = static_cast<int>(invdy * (py - gather_bounds.min[1]));

  if (0 <= x && x < binX && 0 <= y && y < binY) {
    // Add data entries to consecutive bins.
    for (int i = 0; i < min(v.size(), data_width); ++i) {
      addToBin(x, y, i, v[i], false);
    }
    // Increment the bin once.
    incrementBin(x, y);
  }
}

void VolumePlotObject2D::addToBin(RealType px, RealType py, RealType v) {
  addToBin(px, py, 0, v, true);
}

void VolumePlotObject2D::incrementBin(int x, int y) {
  ++counts[x][y];
}
