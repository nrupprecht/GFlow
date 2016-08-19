#include "Theory.h"
#include "ArgParse.h"

int main(int argc, char** argv) {
  // Parameters
  int iters = 100;
  int size = 50;
  double sigma = 0.05;
  double phi = 0.5;
  
  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser(argc, argv);
  pair<string,string> opt;
  stringstream stream;
  opt = parser.find("iters");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> iters;
  }
  opt =parser.find("size");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> size;
  }
  opt=parser.find("points");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> size;
  }
  opt =parser.find("sigma");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> sigma;
  }
  opt =parser.find("phi");
  if (!opt.first.empty()) {
    stream.clear();
    stream << opt.second;
    stream >> phi;
  }
  //----------------------------------------
  
  Theory solver(size, sigma, phi);

  /*
  solver.test("0, 0.952278, 0.000953704, 0.00252778, 0.0543796, 0.606648, 0.14037, \
0.048787, 0.252611, 0.282139, 0.146935, 0.134546, 0.216519, 0.210306, \
0.161398, 0.180824, 0.181935, 0.176315, 0.185361, 0.195157, 0.168102, \
0.173583, 0.188269, 0.198028, 0.19062, 0.180333, 0.184648, 0.191481, \
0.198139, 0.173907, 0.18862, 0.199287, 0.200481, 0.19062, 0.205259, \
0.215389, 0.198139, 0.208278, 0.21687, 0.208852, 0.202472, 0.235037, \
0.223602, 0.22638, 0.217593, 0.239269, 0.249426, 0.246407, 0.253463, \
0.245769, 0.250907, 0.263713, 0.24638, 0.256509, 0.28025, 0.276565, \
0.240593, 0.285667, 0.297148, 0.252037, 0.287796, 0.318898, 0.275917, \
0.27637, 0.338398, 0.284046, 0.277185, 0.363796, 0.294907, 0.286481, \
0.352315, 0.300019, 0.280741, 0.346704, 0.296454, 0.296917, 0.310528, \
0.319676, 0.278, 0.286065, 0.315028, 0.280778, 0.27462, 0.296083, \
0.286102, 0.261944, 0.266861, 0.283324, 0.259222, 0.238852, 0.259917, \
0.269556, 0.211463, 0.220028, 0.250537, 0.227861, 0.205944, 0.221278, \
0.220778, 0.183759, 0.19662, 0.201065, 0.174102, 0.182685, 0.174111, \
0.174731, 0.183407, 0.162963, 0.154037, 0.158491, 0.169685, 0.150981, \
0.157417, 0.150694, 0.142352, 0.141454, 0.139222, 0.149093, 0.131769, \
0.140787, 0.143824, 0.152898, 0.145426, 0.137028, 0.143935, 0.160019, \
0.138231, 0.14938, 0.209398, 0.167269, 0.0888704, 0.166231, 0.435278, \
0.0272963, 0.000185185, 0.000944444, 0.944889, 0", 0.05, 300);
  */

  solver.solve(iters);
  cout << "Time: " << solver.getRunTime() << ";\n\n";

  // Print out all information
  cout << "Print[\"Distribution\"]\n";
  cout << "Dist=" << solver.print() << ";\nListLinePlot[Dist,PlotStyle->Black,PlotRange->All,ImageSize->Large]\n";
  cout << "Print[\"Free Length\"]\n";
  cout << "FL=" << solver.printFreeLength() << ";\nListLinePlot[FL,PlotStyle->Black,ImageSize->Large,PlotRange->{0,1}]\n";
  cout << "Print[\"Collision frequency\"]\n";
  cout << "Freq=" << solver.printFrequency() << ";\nMatrixPlot[Freq,ImageSize->Large]\n";
  /*
  cout << "Print[\"Freq L\"]\n";
  cout << "FreqL=" << solver.printFreqL() << ";\nListLinePlot[FreqL,ImageSize->Large,PlotStyle->Black]\n";
  cout << "Print[\"Freq R\"]\n";
  cout << "FreqR=" << solver.printFreqR() << ";\nListLinePlot[FreqR,ImageSize->Large,PlotStyle->Black]\n";
  */
  cout << "Print[\"Deposition\"]\n";
  cout << "Dep=" << solver.printPropagation() << ";\nMatrixPlot[Dep,ImageSize->Large]\n";
  cout << "Print[\"Difference\"]\n";
  cout << "Diff=" << solver.printDArray() << ";\nListLinePlot[Diff,ImageSize->Large,PlotStyle->Black,PlotRange->All]\n";
  

  solver.printPropagation();

  return 0;
}
