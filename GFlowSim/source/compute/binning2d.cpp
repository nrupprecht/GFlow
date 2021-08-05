#include <compute/binning2d.hpp>

using namespace GFlowSimulation;

Binning2d::Binning2d(int x, int y)
    : dimx(x), dimy(y) {
  if (x <= 0 || y <= 0) {
    throw BadDimension();
  }
  // Allocate bins
  bins = new vector<int>[x * y];
}

Binning2d::~Binning2d() {
  if (bins) {
    delete[] bins;
  }
}

void Binning2d::setBounds(const Bounds &bounds) {
  setBounds(bounds.min[0], bounds.max[0], bounds.min[1], bounds.max[1]);
}

void Binning2d::setBounds(float minx, float maxx, float miny, float maxy) {
  // Set bounds
  min[0] = minx;
  min[1] = miny;
  max[0] = maxx;
  max[1] = maxy;
  // Set bin dimensions
  dx = (maxx - minx) / dimx;
  dy = (maxy - miny) / dimy;
}

void Binning2d::bin_data(const vector<float> &data, int entries, int data_width) {
  // Make sure bounds are good.
  if (min[0] >= max[0] || min[1] >= max[1]) {
    return;
  }
  // Bin data
  for (int i = 0; i < entries; ++i) {
    int bx = (data[i * data_width + 0] - min[0]) / dx;
    int by = (data[i * data_width + 1] - min[1]) / dy;
    if (0 <= bx && bx < dimx && 0 <= by && by < dimy) {
      bins[by * dimx + bx].push_back(i);
    }
  }
}

const vector<int> &Binning2d::at(int x, int y) {
  return bins[y * dimx + x];
}
