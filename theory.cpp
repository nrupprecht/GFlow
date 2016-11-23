#include "ArgParse.h"

#include "NDSolver.h"

double gaussian (gVector vec2d) { return exp(-sqr(3*vec2d.at(0))-sqr(3*vec2d.at(1))); }
double mgaussian(gVector vec2d) { return exp(-sqr(3*vec2d.at(0))-sqr(3*vec2d.at(1)-1)); }
double xgaussian (gVector vec2d) { return exp(-sqr(3*vec2d.at(0))); }

int main(int argc, char** argv) {
  double time = 1;
  double epsilon = 0.01;
  double gridSpacing = 0.05;
  
  //----------------------------------------
  // Parse command line arguments
  //----------------------------------------
  ArgParse parser(argc, argv);
  parser.get("time", time);
  parser.get("epsilon", epsilon);
  parser.get("spacing", gridSpacing);
  //----------------------------------------
  
  NDSolver ndsolver(4);
  ndsolver.setGridSpacing(0, gridSpacing); // X
  ndsolver.setGridSpacing(1, gridSpacing); // Y
  ndsolver.setGridSpacing(2, gridSpacing); // VX
  ndsolver.setGridSpacing(3, gridSpacing); // VY
  ndsolver.setWrapping(0, true);
  ndsolver.setWrapping(1, true);
  ndsolver.setEpsilon(epsilon);
  ndsolver.setField(mgaussian);

  ndsolver.run(time);
  cout << "field=" << mmPreproc(ndsolver.printField()) << ";\n";
  cout << "dfield=" << mmPreproc(ndsolver.printDField()) << ";\n";
  cout << "it=" << ndsolver.getIters() << ";\ntime=" << ndsolver.getTime() << ";\n";
  cout << "ListPlot3D[field,PlotRange->All]\nListPlot3D[dfield,PlotRange->All]";
  return 0;
}
