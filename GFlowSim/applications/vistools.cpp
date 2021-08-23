#include "visualization/visualization.hpp"
#include "compute/field_properties.hpp"

#include "utility/ArgParse.hpp"

using namespace GFlowSimulation;

int main(int argc, char** argv) {
  // --- Seed random numbers
  srand48(time(0));

  // --- Options
  std::string directory = "RunData";
  std::string subdirectory = "general";
  std::string dataName = "Pos-1";
  std::string saveDirectory = "RunData";
  std::string selectionName = "V";
  double radius_multiple = 1.;
  int colorSelect = -1;
  int colorOption = 0;
  int resolution = 1.5*1024;
  // Modes of operation
  bool snapshot = false;
  bool field = false;

  // --- Argument parsing
  ArgParse parser(argc, argv);
  parser.get("directory", directory);
  parser.get("directory", saveDirectory); // By default, save to the same directory
  parser.get("saveDirectory", saveDirectory); // But the defualt can be overruled
  parser.get("selectionName", selectionName);
  parser.get("subdirectory", subdirectory);
  parser.get("data", dataName);
  parser.get("radius_multiple", radius_multiple);
  parser.get("colorSelect", colorSelect);
  parser.get("colorOption", colorOption);
  parser.get("resolution", resolution);
  parser.get("snapshot", snapshot);
  parser.get("field", field);
  // Done finding arguments
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken& illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }

  // --- Create visualization object
  Visualization visualization;
  visualization.setRadiusMultiple(radius_multiple);
  visualization.setColorOption(colorOption);
  visualization.setResolution(resolution);
  if (0<=colorSelect) {
    visualization.setColorSelectionMethod(colorSelect);
    visualization.setSelectionName(selectionName);
  }

  // --- Load the data and create an image
  if (snapshot) {
    visualization.load_and_create(directory+"/general/Snapshot/data.csv", saveDirectory+"/general/Snapshot");
    // Compute field properties
    if (field) {
      FieldProperties field_properties;
      field_properties.load_and_create(directory+"/general/Snapshot/");
    }
  }
  else {
    visualization.load_and_create(directory+"/"+subdirectory+"/"+dataName+"/"+dataName+".csv", saveDirectory+"/"+subdirectory+"/"+dataName);
    // Compute field properties
    if (field) {
      FieldProperties field_properties;
      field_properties.load_and_create(directory+"/"+subdirectory+"/"+dataName);
    }
  }

  // --- End
  return 0;
}
