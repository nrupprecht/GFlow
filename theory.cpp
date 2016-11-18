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

  double time = 1;
  double epsilon = 0.0025;
  double gridSpacing = 0.025;
  
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
  parser.get("epsilon", epsilon);
  parser.get("spacing", gridSpacing);
  //----------------------------------------

  if (bins%2==0) bins++;
  if (vbins%2==0) vbins++; // Make sure we have an odd number of vbins

  NDSolver ndsolver;
  ndsolver.setGridSpacingX(gridSpacing);
  ndsolver.setGridSpacingVX(gridSpacing);
  ndsolver.setEpsilon(epsilon);
  ndsolver.setField(gaussian);

  ndsolver.run(time);
  cout << "field=" << mmPreproc(ndsolver.printField()) << ";\n";
  cout << "dfield=" << mmPreproc(ndsolver.printDField()) << ";\n";
  cout << "it=" << ndsolver.getIters() << ";\ntime=" << ndsolver.getTime() << ";\n";
  cout << "ListPlot3D[field,PlotRange->All]\nListPlot3D[dfield,PlotRange->All]";
  return 0;
}
