#include "persistence-length.hpp"

namespace GFlowSimulation {

  PersistenceLength::PersistenceLength(GFlow *gflow) : GraphObject(gflow, "PersistenceLength", "N links away", "<R0*Rk>"), Group(gflow) {};

  PersistenceLength::PersistenceLength(GFlow *gflow, Group& group) : GraphObject(gflow, "PersistenceLength", "N links away", "<R0*Rk>"), Group(group) {};

  void PersistenceLength::pre_integrate() {
    // Set up bins.
    makeBins(1., Group::size()+1, Group::size()-2);
    // Reset.
    ndata_points = 0;
  }

  void PersistenceLength::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check() || Group::size()<3) return;

    // Update local ids?
    if (locals_changed) update_local_ids();
    locals_changed = false;
   
    // Get the position array. Create helper vectors.
    auto x = simData->X();
    Vec dX1(sim_dimensions), dX2(sim_dimensions);
    RealType projection = 0;

    // First angle.
    int id1 = at(1), id0 = at(0);
    gflow->getDisplacement(x(id1), x(id0), dX1.data);
    dX1.normalize();

    // Calculate the average R0*Rk for bonds.
    for (int i=1; i<Group::size()-1; ++i) {
      id0 = at(i);
      id1 = at(i+1);

      // Find the bonds.
      gflow->getDisplacement(x(id1), x(id0), dX2.data);
      dX2.normalize();
      // Bin data.
      data[i-1].second += dX2*dX1;
    }
    // Increment number of data points.
    ++ndata_points;
  }

  void PersistenceLength::post_integrate() {
    if (ndata_points>0) 
      for (int i=0; i<data.size(); ++i) {
        data[i].first = i+1;
        data[i].second /= ndata_points;
      }
  }

}