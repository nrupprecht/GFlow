#include "ArgParse.h"

#include "NDSolver.h"

double gaussian (gVector vec2d) { return exp(-sqr(3*vec2d.at(0))-sqr(3*vec2d.at(1))); }
double mgaussian(gVector vec2d) { return exp(-sqr(3*vec2d.at(0))-sqr(3*vec2d.at(1)-1)); }
double xgaussian (gVector vec2d) { return exp(-sqr(3*vec2d.at(0))); }

int main(int argc, char** argv) {
  // Start the clock
  auto start_t = clock();

  // Parameters
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
  // Set X
  ndsolver.setBounds(0, -1, 1, false);
  ndsolver.setGridSpacing(0, gridSpacing, false);
  // Set Y
  ndsolver.setBounds(1, -1, 1, false);
  ndsolver.setGridSpacing(1, gridSpacing, false);
  // Set VX 
  ndsolver.setBounds(2, -1, 1, false);
  ndsolver.setGridSpacing(2, gridSpacing, false);
  // Set VY
  ndsolver.setBounds(3, -1, 1, false);
  ndsolver.setGridSpacing(3, gridSpacing, false); // VY
  // Remake fields
  ndsolver.remake();
  
  ndsolver.setWrapping(0, true);
  ndsolver.setWrapping(1, true);
  ndsolver.setEpsilon(epsilon);
  ndsolver.setField(mgaussian);

  cout << "----------------------- RUN SUMMARY -----------------------\n\n";
  cout << "Command: ";
  for (int i=0; i<argc; i++) cout << argv[i] << " ";


  cout << "\n...........................................................\n\n";

  ndsolver.run(time);
  auto end_t = clock(); // End timing

  /// Print the run information
  cout << "Sim Time: " << time << ", Run time: " << ndsolver.getRunTime() << ", Ratio:\
 " << time/ndsolver.getRunTime() << endl;
  cout << "Actual (total) program run time: " << (double)(end_t-start_t)/CLOCKS_PER_SEC \
       << "\n";
  cout << "Iterations: " << ndsolver.getIters();
  cout << "\n\n----------------------- END SUMMARY -----------------------\n\n";

  cout << "field=" << mmPreproc(ndsolver.printField()) << ";\n";
  cout << "dfield=" << mmPreproc(ndsolver.printDField()) << ";\n";
  cout << "it=" << ndsolver.getIters() << ";\ntime=" << ndsolver.getTime() << ";\n";
  cout << "ListPlot3D[field,PlotRange->All]\nListPlot3D[dfield,PlotRange->All]";
  return 0;
}
