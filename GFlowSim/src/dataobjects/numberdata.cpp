#include "numberdata.hpp"
// Other files
#include "../visualization/palette.hpp"

namespace GFlowSimulation {

  NumberData::NumberData(GFlow *gflow) : DataObject(gflow, "NumberData") {};

  void NumberData::post_step() {
    int ntypes = simData->ntypes();

    if (numberData.empty()) numberData = vector<vector<RPair> >(ntypes);
    vector<int> nums(ntypes, 0);

    int size = simData->size_owned();
    for (int n=0; n<size; ++n) {
      int type = simData->Type(n);
      if (-1<type && type<ntypes)
        ++nums[type];
    }

    // Store data
    RealType time = gflow->getElapsedTime();
    for (int i=0; i<ntypes; ++i)
      numberData[i].push_back(RPair(time, static_cast<RealType>(nums[i])));
  }

  bool NumberData::writeToFile(string fileName, bool) {

    // Draw a graph using a palette object
    Palette graph(1024,512);
    GraphOptions options;
    options.setMinY(0);
    options.setMaxY(1500);
    options.useBcgd = true;

    for (auto &data : numberData) {
      RGBApixel line = randomColor();
      options.setLineColor(line);
      graph.drawGraph2d(data, options);
    }
    graph.writeToFile(fileName+"/numberdata.bmp");

    return true;
  }


}