#include "groupnetforce.hpp"

namespace GFlowSimulation {

  GroupNetForce::GroupNetForce(GFlow *gflow) : MultiGraphObject(gflow, "GroupForce", "time", "force", gflow->getSimDimensions()), Group(gflow) 
  {
    for (int i=0; i<sim_dimensions; ++i)
      axis_y[i] = "Force - F[" + toStr(i) + "]";
  };

  void GroupNetForce::post_step() {
    // Only record if enough time has gone by and there are particles to keep track of.
    if (!DataObject::_check()) return;

    // Make sure local ids are up to date.
    if (locals_changed) {
      update_local_ids();
      locals_changed = false;
    }

    // Add a new entry to modify.
    addEntry();
    // Set the time.
    //getX() = gflow->getElapsedTime();
    // Set the forces.
    Vec F(sim_dimensions);
    findNetForce(F.data);

    // Gather and store data on processor 0.
    gatherData(gflow->getElapsedTime(), F);
    /*
    // Data reduction.
    MPIObject::mpi_sum0(F.data, sim_dimensions);

    // Set each dimension of force.
    if (topology->getRank()==0) {
      for (int d=0; d<sim_dimensions; ++d) getY(d) = F[d];
    }
    */
  }

}
