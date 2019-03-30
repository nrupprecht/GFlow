#include "multigraphobject.hpp"

namespace GFlowSimulation {

  MultiGraphObject::MultiGraphObject(GFlow *gflow, const string& name, int nd) : DataObject(gflow, name, DataObjectType::MULTIGRAPH), ndata(nd) {
    // Must have ndata >= 1
    if (ndata<=0) throw EmptyData("Must have ndata>=1");
    // Set up data array.
    resetData();
  };

  MultiGraphObject::MultiGraphObject(GFlow *gflow, const string& name, const string& ax, const string& ay, int nd) 
    : DataObject(gflow, name, DataObjectType::MULTIGRAPH), axis_x(ax), axis_y(ay), ndata(nd) 
  {
    // Must have ndata >= 1
    if (ndata<=0) throw EmptyData("Must have ndata>=1");
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
    ofstream fout(dirName+dataName+".csv");
    if (fout.fail()) return false;
    // Print out the data
    for (int i=0; i<ndata_points; ++i) {
      for (int j=0; j<ndata+1; ++j) {
        fout << multi_data[j][i];
        if (j!=ndata) fout << ",";
      }
      fout << endl;
    }

    fout.close();
    // Print out the axes labels
    fout.open(dirName+"axes.csv");
    if (fout.fail()) return false;
    fout << axis_x << endl << axis_y << endl;
    fout.close();
    
    /*
    // Optionally print a plot of the data using vistools.
    if (print_plot) {
      // Draw a graph using a palette object
      Palette graph(1024,512);
      GraphOptions options;
      options.setMinY(0);
      options.setBackground(RGB_White);
      options.setLineColor(RGB_Red);
      // Draw the graph
      graph.drawGraph2d(data, options);
      graph.writeToFile(dirName+"/"+dataName+".bmp");
    }
    */

    // Return success
    return true;
  }

  void MultiGraphObject::resetData() {
    // Clear the data vector by recreating it.
    multi_data = vector< vector<RealType> >(ndata+1, vector<RealType>());
    // Reset number of data points counter.
    ndata_points = 0;
  }

  void MultiGraphObject::addEntry() {
    for (int i=0; i<multi_data.size(); ++i) {
      multi_data[i].push_back(0);
    }
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

}