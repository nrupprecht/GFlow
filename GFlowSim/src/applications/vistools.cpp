#include "../visualization/visualization.hpp"

#include "../utility/ArgParse.hpp"

using namespace GFlowSimulation;

int main(int argc, char** argv) {

  // --- Seed random numbers
  srand48(time(0));

  // --- Options

  string directory = "RunData";
  string subdirectory = "Pos";
  string saveDirectory = "RunData";
  double radius_multiple = 1.;
  int colorOption = 0;
  int resolution = 1.5*1024;

  // Argument parsing
  ArgParse parser(argc, argv);
  parser.get("directory", directory);
  parser.get("directory", saveDirectory); // By default, save to the same directory
  parser.get("saveDirectory", saveDirectory); // But the defualt can be overruled
  parser.get("subdirectory", subdirectory);
  parser.get("radius_multiple", radius_multiple);
  parser.get("colorOption", colorOption);
  parser.get("resolution", resolution);
  // Done finding arguments
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }

  // Create visualization object
  Visualization visualization;
  visualization.setRadiusMultiple(radius_multiple);
  visualization.setColorOption(colorOption);
  visualization.setResolution(resolution);

  visualization.load_and_create(directory+"/"+subdirectory+"/data.csv", saveDirectory+"/"+subdirectory);


  return 0;
}
