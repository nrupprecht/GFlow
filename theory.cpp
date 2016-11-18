#include "Checker.h"
#include "ArgParse.h"

#include "NDSolver.h"

double gaussian (gVector vec2d) { return exp(-sqr(3*vec2d.at(0))-sqr(3*vec2d.at(1))); }
double xgaussian (gVector vec2d) { return exp(-sqr(3*vec2d.at(0))); }

int main(int argc, char** argv) {
  // Parameters 
  int iters = 100;
  int bins = 31;
  int vbins = 21;
  int recDelay = 10;
  double radius = 0.05;
  string filename= "";

  double time = 0.05;
  
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

  parser.get("time", time);
  //----------------------------------------

  if (bins%2==0) bins++;
  if (vbins%2==0) vbins++; // Make sure we have an odd number of vbins

  NDSolver ndsolver;
  ndsolver.setGridSpacingX(0.05);
  ndsolver.setGridSpacingVX(0.05);
  ndsolver.setEpsilon(0.05);
  ndsolver.setField(gaussian);
  //ndsolver.setField(xgaussian);
  ndsolver.run(time);
  cout << "field=" << mmPreproc(ndsolver.printField()) << ";\n";
  cout << "dfield=" << mmPreproc(ndsolver.printDField()) << ";\n";
  cout << "it=" << ndsolver.getIters() << ";\ntime=" << ndsolver.getTime() << ";\n";
  cout << "ListPlot3D[field,PlotRange->All]\nListPlot3D[dfield,PlotRange->All]";
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
