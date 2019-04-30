#include "graphobject.hpp"
// Other files
#include "../../visualization/palette.hpp"

namespace GFlowSimulation {

  GraphObject::GraphObject(GFlow *gflow, const string& name) : DataObject(gflow, name, DataObjectType::GRAPH) {};

  GraphObject::GraphObject(GFlow *gflow, const string& name, const string& ax, const string& ay) 
    : DataObject(gflow, name, DataObjectType::GRAPH), axis_x(ax), axis_y(ay) {};

  void GraphObject::pre_integrate() {
    // Parent class.
    DataObject::pre_integrate();
    // Clear data.
    data.clear();
  }

  void GraphObject::addEntry(RealType t, RealType f) {
    data.push_back(RPair(t, f));
  }

  RealType GraphObject::ave() const {
    if (data.empty()) return 0;
    // Accumulator
    RealType total = 0;
    // Tally 
    for (auto v : data) total += v.second;
    // Temporal average
    return total / (data.size());
  }

  void GraphObject::setAxes(const string& ax, const string& ay) {
    axis_x = ax;
    axis_y = ay;
  }

  bool GraphObject::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (data.empty()) return true;
    // The name of the directory for this data
    string dirName = _correctDirName(fileName);

    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);
    ofstream fout(dirName+dataName+"-"+toStr(object_counter)+".csv");
    if (fout.fail()) return false;
    // Print out the data
    for (auto v : data) fout << v.first << "," << v.second << endl;
    fout.close();
    // Print out the axes labels
    fout.open(dirName+"axes.csv");
    if (fout.fail()) return false;
    fout << axis_x << endl << axis_y << endl;
    fout.close();
    
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

    // Return success
    return true;
  }

}