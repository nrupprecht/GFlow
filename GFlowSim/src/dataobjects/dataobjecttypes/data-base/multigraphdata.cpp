#include "multigraphdata.hpp"

namespace GFlowSimulation {

  MultiGraphData::MultiGraphData(int nd) : ndata(nd) {
    // Must have ndata >= 1
    if (ndata<=0) throw EmptyData("Must have ndata>=1");
    // Write data flags are true by default.
    write_data = vector<bool>(nd, true);
    // Default y axis name.
    axis_y = vector<string>(nd, "y");
    // Set up data array.
    resetData();
  };

  MultiGraphData::MultiGraphData(const string& ax, const string& ay, int nd) : axis_x(ax), ndata(nd) {
    // Must have ndata >= 1
    if (ndata<=0) throw EmptyData("Must have ndata>=1");
    // Write data flags are true by default.
    write_data = vector<bool>(nd, true);
    // All y axes the same.
    axis_y = vector<string>(nd, ay);
    // Set up data array.
    resetData();
  };

  bool MultiGraphData::write_to_file(const string& dirName, const string& dataName) {
    // Check if there's anything to do
    if (multi_data.empty() || multi_data[0].empty()) return true;

    ofstream fout(dirName+dataName+".csv");
    if (fout.fail()) {
      cout << "File [" << (dirName+"/"+dataName+".csv") << "] failed to open for MultiGraphData write.\n";
      return false;
    }
    // Print out the data
    for (int i=0; i<multi_data[0].size(); ++i) {
      for (int j=0; j<ndata+1; ++j) {
        if (j==0 || write_data[j-1]) {
          // Prepend with a "," (if not the first entry) so there will not be null entries if the last entry should not be printed.
          if (j!=0) fout << ",";
          fout << multi_data[j][i];
        }
      }
      fout << endl;
    }
    fout.close();
    // Print out all the axes labels.
    fout.open(dirName+"/axes.csv");
    if (fout.fail()) {
      cout << "File [" << (dirName+"/axes.csv") << "] failed to open for MultiGraphData write.\n";
      return false;
    }
    for (int i=0; i<ndata; ++i)
      if (write_data[i]) fout << axis_x << "," << axis_y[i] << endl;
    fout.close();

    // Return success
    return true;
  }

  real MultiGraphData::ave(int i) {
    if (i<0 || ndata<=i || multi_data.empty() || multi_data[0].size()==0) return 0;
    real total = 0;
    // Tally 
    for (auto v : multi_data[i+1]) total += v;
    // Return the average
    return total/multi_data[0].size();
  }

  vector<RPair> MultiGraphData::getEntry(int i) const {
    if (i<0 || ndata<=i || multi_data.empty() || multi_data[0].size()==0) return vector<RPair>();
    // Accumulate data
    vector<RPair> entry;
    for (int j=0; j<multi_data[0].size(); ++j)
      entry.push_back(RPair(multi_data[0][j], multi_data[i+1][j]));
    // Return data
    return entry;
  }

  int MultiGraphData::size() const {
    return multi_data.size();
  }

  void MultiGraphData::resetData(int size) {
    if (0<size) multi_data = vector< vector<real> >(ndata+1, vector<real>(size, 0));
    else multi_data = vector< vector<real> >(ndata+1, vector<real>());
  }

  void MultiGraphData::addEntry() {
    for (int i=0; i<multi_data.size(); ++i) {
      multi_data[i].push_back(0);
    }
  }

  void MultiGraphData::addEntry(real v) {
    for (int i=0; i<multi_data.size(); ++i) {
      multi_data[i].push_back(0);
    }
    // Set X value.
    multi_data[0][multi_data[0].size()-1] = v;
  }

  real& MultiGraphData::getX() {
    // Check for multi data being empty.
    if (multi_data[0].empty()) throw EmptyData();
    // Otherwise, there is data.
    return multi_data[0][multi_data[0].size()-1];
  }

  real& MultiGraphData::getY(int i) {
    // Check for multi data being empty.
    if (multi_data[0].empty()) throw EmptyData();
    // Otherwise, there is data. Check if we are in bounds.
    if (i<0 || ndata<=i) throw MultiDataOutOfBounds("Out of bounds: i="+toStr(i));
    // Otherwise, return data.
    return multi_data[i+1][multi_data[i+1].size()-1];
  }

  real& MultiGraphData::atX(int i) {
    return multi_data[0][i];
  }

  real& MultiGraphData::atY(int d, int i) {
    return multi_data[1+d][i];
  }

  void MultiGraphData::gatherAverageData(const real x, const Vec yvals, int count) {
    // Gather data and data counts on processor 0.
    MPIObject::mpi_sum0(count);
    MPIObject::mpi_sum0(yvals.data, yvals.size());
    // If this is the root processor (0), store the data.
    if (MPIObject::getRank()==0) {
      // Add an entry with the given x value.
      addEntry(x);
      // Store data
      real c = static_cast<real>(count);
      for (int d=0; d<yvals.size(); ++d) getY(d) = yvals[d] / c;
    }
  }

  void MultiGraphData::gatherData(const real x, const Vec yvals) {
    // Gather data and data counts on processor 0.
    MPIObject::mpi_sum0(yvals.data, yvals.size());
    // If this is the root processor (0), store the data.
    if (MPIObject::getRank()==0) {
      // Add an entry with the given x value.
      addEntry(x);
      // Store data
      for (int d=0; d<yvals.size(); ++d) getY(d) = yvals[d];
    }
  }

}