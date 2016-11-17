#include "Checker.h"
#include "ArgParse.h"

#include "NDSolver.h"

int main(int argc, char** argv) {
  // Parameters 
  int iters = 100;
  int bins = 31;
  int vbins = 21;
  int recDelay = 10;
  double radius = 0.05;
  string filename= "";
  
  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser(argc, argv);
  parser.get("iters", iters);
  parser.get("bins", bins);
  parser.get("vbins", vbins);
  parser.get("recDelay", recDelay);
  parser.get("radius", radius);
  parser.get("filename", filename);
  //----------------------------------------

  if (bins%2==0) bins++;
  if (vbins%2==0) vbins++; // Make sure we have an odd number of vbins


  NDSolver ndsolver;
  ndsolver.run(1.);
  return 0;

  Checker checker(bins, vbins, radius);
  checker.setInitialCondition();
  checker.setFlowV(0.2);
  checker.setEpsilon(0.001);
  checker.setRecDelay(recDelay);
  if (filename!="") 
    if (!checker.readFromFile(filename)) {
      cout << "Failed to read from file.";
      return 0;
    }

  // Run
  checker.run(iters);

  stringstream stream;
  string str;

  // Profile
  stream << checker.getProfile();
  stream >> str;
  cout << "profile=" << mmPreproc(str) << ";\n";
  cout << "Print[\"Profile\"]\n";
  cout << "ListLinePlot[profile,ImageSize->Large,PlotRange->All]\n";
  stream.clear();
  str.clear();

  stream << checker.getDProfileDT();
  stream>> str;
  cout << "dprofile=" << mmPreproc(str) << ";\n";
  cout << "Print[\"Time derivative of profile\"]\n";
  cout << "ListLinePlot[dprofile,ImageSize->Large,PlotRange->All]\n";
  stream.clear();
  str.clear();

  stream << checker.getSpaceProfile();
  stream >> str;
  cout << "fullprofile=" << mmPreproc(str) << ";\n";
  cout << "Print[\"Full Profile\"]\n";
  cout << "ArrayPlot[fullprofile,ImageSize->Large]\n";
  stream.clear();
  str.clear();

  stream << checker.getSpaceAnimation();
  stream >> str;
  cout << "animation=" << mmPreproc(str) << ";\n";
  stream.clear();
  str.clear();
  cout << "Print[\"Animation\"]\n";
  cout << "ListAnimate[Table[ArrayPlot[animation[[i]],ImageSize->Large],{i,1,Length[animation]}]]\n";

  stream << checker.getCollapsedDistribution();
  stream >> str;
  cout << "dist=" << mmPreproc(str) << ";\n";
  cout << "Print[\"Collapsed Distribution (y is collapsed, cycles through x)\"]\n";
  cout << "ListAnimate[Table[MatrixPlot[dist[[i]]],{i,1,Length[dist]}]]\n";
  stream.clear();
  str.clear();

  return 0;
}
