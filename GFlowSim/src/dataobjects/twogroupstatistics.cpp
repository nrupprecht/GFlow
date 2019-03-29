#include "twogroupstatistics.hpp"

namespace GFlowSimulation {

  TwoGroupStatistics::TwoGroupStatistics(GFlow *gflow) : DataObject(gflow, "TwoGroupStatistics") {
    setBins(50, 101);
  };

  void TwoGroupStatistics::pre_integrate() {
    setBins(bins_x, bins_y);
  }

  void TwoGroupStatistics::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check() || group1.size()==0 || group2.size()==0) return;

    // Find the local ids of the group atoms if simdata has altered the local ids.
    if (locals_changed) {
      group1.update_local_ids(simData);
      group2.update_local_ids(simData);
    }
    locals_changed = false;

    // Find the centers of mass for each group
    RealType *rcm1 = new RealType[sim_dimensions], *rcm2 = new RealType[sim_dimensions], *dr = new RealType[sim_dimensions];
    group1.findCenterOfMass(rcm1, simData);
    group2.findCenterOfMass(rcm2, simData);
    gflow->getDisplacement(rcm2, rcm1, dr);
    RealType R = sqrt(dotVec(dr, dr, sim_dimensions));

    // If within the cutoff radius, bin the data
    if (R<r_max) {
      // Normalize dr
      normalizeVec(dr, sim_dimensions);
      // Find the net force on each group
      RealType *f1 = new RealType[sim_dimensions], *f2 = new RealType[sim_dimensions], *dF = new RealType[sim_dimensions];
      group1.findNetForce(f1, simData);
      group2.findNetForce(f2, simData);
      subtractVec(f2, f1, dF, sim_dimensions);
      RealType rel_f = dotVec(dF, dr, sim_dimensions);
      // Find the com velocity for each group
      RealType *v1 = new RealType[sim_dimensions], *v2 = new RealType[sim_dimensions], *vr = new RealType[sim_dimensions];
      group1.findCOMVelocity(v1, simData);
      group2.findCOMVelocity(v2, simData);
      subtractVec(v2, v1, vr, sim_dimensions);
      RealType rel_v = dotVec(vr, dr, sim_dimensions);
      
      // Distance bin
      int br = (R/r_max)*bins_x;
      // Force bin
      int bf = static_cast<int>(0.5*bins_y*(rel_f + f_max)/f_max);
      if (0<=bf && bf<bins_y) ++count_f[br][bf];
      
      // Velocity bin
      int bv = static_cast<int>(0.5*bins_y*(rel_v + v_max)/v_max);
      if (0<=bv && bv<bins_y) ++count_v[br][bv];

      // Clean up
      delete [] f1;
      delete [] f2;
      delete [] dF;
      delete [] v1;
      delete [] v2;
      delete [] vr;
    }

    // Clean up
    delete [] rcm1;
    delete [] rcm2; 
    delete [] dr;
   
  }

  bool TwoGroupStatistics::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (count_v.empty() && count_f.empty()) return true;
    // The name of the directory for this data
    string dirName = _correctDirName(fileName);

    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);
    // Open file for velocity data
    ofstream fout(dirName+dataName+"-vel.csv");
    if (fout.fail()) return false;
    // Print out the velocity data
    RealType dr = r_max/bins_x;
    RealType dv = 2*v_max/bins_y;
    RealType r = 0.5*dr;
    for (int x=0; x<bins_x; ++x) {
      RealType v = -v_max;
      for (int y=0; y<bins_y; ++y) {
        fout << r << "," << v << "," << count_v[x][y] << endl;
        v += dv;
      }
      r += dr;
    }
    fout.close();

    // Open file for force data
    fout.open(dirName+dataName+"-force.csv");
    if (fout.fail()) return false;
    // Print out the velocity data
    RealType df = 2*f_max/(bins_y-1);
    r = 0;
    for (int x=0; x<bins_x; ++x) {
      RealType f = -f_max;
      for (int y=0; y<bins_y; ++y) {
        fout << r << "," << f << "," << count_f[x][y] << endl;
        f += df;
      }
      r += dr;
    }
    fout.close();

    // Return success
    return true;
  }

  void TwoGroupStatistics::setBins(int x, int y) {
    // Set values
    bins_x = x;
    bins_y = y;
    // Set arrays
    count_v = vector<vector<int> >(x, vector<int>(y, 0));
    count_f = vector<vector<int> >(x, vector<int>(y, 0));
  }

  void TwoGroupStatistics::setGroups(const Group& g1, const Group& g2) {
    group1 = g1;
    group2 = g2;
  }

}