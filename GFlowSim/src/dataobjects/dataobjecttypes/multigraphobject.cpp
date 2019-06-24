#include "multigraphobject.hpp"

namespace GFlowSimulation {

  MultiGraphObject::MultiGraphObject(GFlow *gflow, const string& name, int nd) : DataObject(gflow, name, DataObjectType::MULTIGRAPH), ndata(nd) {
    // Must have ndata >= 1
    if (ndata<=0) throw EmptyData("Must have ndata>=1");
    // Write data flags are true by default.
    write_data = vector<bool>(nd, true);
    // Default y axis name.
    axis_y = vector<string>(nd, "y");
    // Set up data array.
    resetData();
  };

  MultiGraphObject::MultiGraphObject(GFlow *gflow, const string& name, const string& ax, const string& ay, int nd) 
    : DataObject(gflow, name, DataObjectType::MULTIGRAPH), axis_x(ax), ndata(nd) 
  {
    // Must have ndata >= 1
    if (ndata<=0) throw EmptyData("Must have ndata>=1");
    // Write data flags are true by default.
    write_data = vector<bool>(nd, true);
    // All y axes the same.
    axis_y = vector<string>(nd, ay);
    // Set up data array.
    resetData();
  };

  void MultiGraphObject::pre_integrate() {
    // Call parent class.
    DataObject::pre_integrate();
    // Reset data
    resetData();
  }

  bool MultiGraphObject::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (multi_data.empty() || ndata_points==0) return true;
    // The name of the directory for this data
    string dirName = _correctDirName(fileName);
    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);
    ofstream fout(dirName+dataName+"-"+toStr(object_counter)+".csv");
    if (fout.fail()) return false;
    // Print out the data
    for (int i=0; i<ndata_points; ++i) {
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
    fout.open(dirName+"axes.csv");
    if (fout.fail()) return false;
    for (int i=0; i<ndata; ++i)
      if (write_data[i]) fout << axis_x << "," << axis_y[i] << endl;
    fout.close();

    // Return success
    return true;
  }

  RealType MultiGraphObject::ave(int i) {
    if (i<0 || ndata<=i || ndata_points==0) return 0;
    RealType total = 0;
    // Tally 
    for (auto v : multi_data[i+1]) total += v;
    // Return the average
    return total/ndata_points;
  }

  vector<RPair> MultiGraphObject::getEntry(int i) {
    if (i<0 || ndata<=i || ndata_points==0) return vector<RPair>();
    // Accumulate data
    vector<RPair> entry;
    for (int j=0; j<ndata_points; ++j)
      entry.push_back(RPair(multi_data[0][j], multi_data[i+1][j]));
    // Return data
    return entry;
  }

  void MultiGraphObject::resetData(int size) {
    // Clear the data vector by recreating it.
    if (size>0) {
      multi_data = vector< vector<RealType> >(ndata+1, vector<RealType>(size, 0));
      ndata_points = size;
    }
    else {
      multi_data = vector< vector<RealType> >(ndata+1, vector<RealType>());
      ndata_points = 0;
    }
  }

  void MultiGraphObject::addEntry() {
    for (int i=0; i<multi_data.size(); ++i) {
      multi_data[i].push_back(0);
    }
    // Increment counter
    ++ndata_points;
  }

  void MultiGraphObject::addEntry(RealType v) {
    for (int i=0; i<multi_data.size(); ++i) {
      multi_data[i].push_back(0);
    }
    // Set X value.
    multi_data[0][ndata_points] = v;
    // Increment counter
    ++ndata_points;
  }

  RealType& MultiGraphObject::getX() {
    // Check for multi data being empty.
    if (multi_data[0].empty()) throw EmptyData();
    // Otherwise, there is data.
    return multi_data[0][ndata_points-1];
  }

  RealType& MultiGraphObject::getY(int i) {
    // Check for multi data being empty.
    if (multi_data[0].empty()) throw EmptyData();
    // Otherwise, there is data. Check if we are in bounds.
    if (i<0 || ndata<=i) throw MultiDataOutOfBounds("Out of bounds: i="+toStr(i));
    // Otherwise, return data.
    return multi_data[i+1][ndata_points-1];
  }

  RealType& MultiGraphObject::atX(int i) {
    return multi_data[0][i];
  }

  RealType& MultiGraphObject::atY(int d, int i) {
    return multi_data[1+d][i];
  }

}