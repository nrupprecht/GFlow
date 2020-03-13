#include "cavity-positions.hpp"

namespace GFlowSimulation {

  CavityPositions::CavityPositions(GFlow *gflow, real v) : PositionData(gflow), MultiGraphData("time", "", 3), limit_velocity(v) {
    dataName = "Cav";
    select_function = [=](shared_ptr<SimData> simdata, int n)->bool {
      return simdata->V(n, 0)<0.5f*limit_velocity && 0<simdata->X(n, 0) && 0<simdata->Im(n);
    };
  };

  void CavityPositions::post_step() {
    // This will gather all the particle data, which we will then use to record multigraph data.
    PositionData::post_step();

    if (topology->getRank()==0) {
      Vec data(3);
      auto &particle_data = positions[positions.size()-1];
      int data_width = storeData.getDataWidth(), num_particles = particle_data.size() / data_width;
      data[0] = num_particles;
      for (int i=0; i<num_particles; ++i) {
        data[1] += particle_data[i*data_width];
        data[2] += particle_data[i*data_width+1];
      }
      data[1] /= num_particles;
      data[2] /= num_particles;
    }
  }

  bool CavityPositions::writeToFile(string fileName, bool make_directory) {
    // The name of the directory for this data
    string dirName = _correctDirName(fileName);
    // Create a directory for all the data
    string name = dirName+dataName+"-"+toStr(object_counter)+".csv";
    return PositionData::writeToFile(fileName) && write_to_file(name);
  }

}