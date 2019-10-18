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

  RPair GraphObject::first() const {
    if (!data.empty()) return data[0];
    else return RPair(0., 0.);
  }
  
  RPair GraphObject::last() const {
    if (!data.empty()) return data[data.size()-1];
    else return RPair(0., 0.);
  }

  int GraphObject::size() const {
    return data.size();
  }

  void GraphObject::setAxes(const string& ax, const string& ay) {
    axis_x = ax;
    axis_y = ay;
  }

  void GraphObject::makeBins(const RealType min, const RealType max, const int nbins) {
    data.clear();
    if (nbins>0) {
      data = vector<RPair>(nbins, RPair(0,0));
      // Set bin ticks.
      if (nbins>1) {
        RealType dbin = (max-min)/(nbins-1);
        for (int i=0; i<nbins; ++i) data[i].first = min + dbin*i;
      }
      else data[0].first = 0.5*(min+max);
    }
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
    
    // Return success
    return true;
  }

  void GraphObject::gatherAverageData(const RealType x, RealType yval, int count) {
    // Gather data and data counts on processor 0.
    MPIObject::mpi_sum0(count);
    MPIObject::mpi_sum0(yval);
    // Check for NANs
    if (isnan(yval)) throw NanValue("GraphObject::gatherAverageData");
    // If this is the root processor (0), store the data average.
    if (topology->getRank()==0) addEntry(x, yval/count);
  }

  void GraphObject::gatherAverageData(const RealType x, RealType yval, RealType weight) {
    // Gather data and data counts on processor 0.
    MPIObject::mpi_sum0(weight);
    MPIObject::mpi_sum0(yval);
    // Check for NANs
    if (isnan(yval)) throw NanValue("GraphObject::gatherAverageData");
    // If this is the root processor (0), store the data average.
    if (topology->getRank()==0) addEntry(x, yval/weight);
  }

  void GraphObject::gatherData(const RealType x, RealType yval) {
    MPIObject::mpi_sum0(yval);
    // Check for NANs
    if (isnan(yval)) throw NanValue("GraphObject::gatherAverageData");
    // If this is the root processor (0), store the data.
    if (topology->getRank()==0) addEntry(x, yval);
  }

}
