#include "demon.hpp"
// Other files
#include "../base/domainbase.hpp"
#include "../alldataobjects.hpp"
#include "../base/datamaster.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  Demon::Demon(GFlow *gflow) : Modifier(gflow), tau(0.1) {
    // Set up data objects
    kineticL = new GraphObject(gflow, "KineticL", "t", "KE");
    kineticR = new GraphObject(gflow, "KineticR", "t", "KE");
    numberL  = new GraphObject(gflow, "NumberL", "t", "number");
    numberR  = new GraphObject(gflow, "NumberR", "t", "number");
    
    // Add data objects to gflow.
    gflow->addDataObject(kineticL);
    gflow->addDataObject(kineticR);
    gflow->addDataObject(numberL);
    gflow->addDataObject(numberR);

    // Create parameters object
    parameters = new Parameters(gflow);
    gflow->addDataObject(parameters);

    // Default door choice function
    check_choice = &direction_demon;
  };

  void Demon::pre_integrate() {
    // Only run in simulation mode
    if (gflow->getRunMode()!=RunMode::SIM) return;

    // Get the side
    side_entry = simData->requestIntegerData("Side");
    // Entry must be positive. Otherwise, something is wrong.
    if (side_entry<0) throw false;
    // Assign particles their side.
    assign_side(false);
    // Check point the initial data.
    checkpoint_data();
    // Set the last check time.
    last_check = gflow->getElapsedTime();

    // --- Record where door particles are

    // Set up vector of door particle positions.
    door_positions = vector<Vec>(size(), Vec(sim_dimensions));
    // Update local ids.
    update_local_ids(simData);
    // Get data from simdata.
    RealType **x = simData->X();
    // Put particles back
    for (int i=0; i<size(); ++i) {
      int id = at(i);
      // Put the particle back in place.
      copyVec(x[id], door_positions[i]);
    }

    // --- Look for an animation object, so we can modify it.
    for (auto obj : dataMaster->getDataObjects()) {
      PositionData *ob = dynamic_cast<PositionData*>(obj);
      if (ob) animate_object = ob;
    }

    // Open the door
    open_door();
  }

  void Demon::pre_forces() {
    // Only run in simulation mode
    if (gflow->getRunMode()!=RunMode::SIM) return;

    // Get the current time
    RealType time = gflow->getElapsedTime();

    // Check if enough time has gone by.
    if (time - last_check >= tau) {
      // If door was closed, then it was previously decided that the last time window should have been a closed window.
      // In this case, nl = nr = el = er = 0, and we can setup to check the next time window.
      // If the door was open, then we should count the changes and decide if we need to go back in time and keep the door
      // closed or not.
      if (door_open) {

        int nl = 0, nr = 0;
        RealType el = 0, er = 0;
        // What happened during the last time window.
        count_changes(nl, nr, el, er);
        // Check whether door should open
        bool open = should_door_open(nl, nr, el, er);
        // If it was good that the door was open, then move on.
        if (open) {
          add_data(nl - nr, el - er);
          // Assign the particle's sides.
          assign_side();
          // Record the door time point.
          last_check = time;
          // Checkpoint
          checkpoint_data();
          // Add data to data objects.
          kineticL->addEntry(time, El);
          kineticR->addEntry(time, Er);
          numberL->addEntry(time, Nl);
          numberR->addEntry(time, Nr);
          
          // The door is still open, so we don't have to open it.
        }
        // Otherwise, go back in time and keep the door closed.
        else {
          // Revert to last checkpoint. The particles' side data will be correct after this.
          revert_to_last_checkpoint(false); // Don't rebuild, since closing the door will rebuild.
          // Close the door.
          close_door();
          // Do not checkpoint data or assign side.
          assign_side();
        }
      }
      // Door was closed.
      else {
        // Door was closed, so nl = nr = el = er = 0. Store this.
        add_data(0, 0);
        // Add data to data objects.
        kineticL->addEntry(time, El);
        kineticR->addEntry(time, Er);
        numberL->addEntry(time, Nl);
        numberR->addEntry(time, Nr);
        // Record the door time point.
        last_check = time;
        // Checkpoint
        checkpoint_data();
        // Open the door to see what will happen.
        open_door();
        // The door was closed, so there is no need to reassign sided.
      }
    }
  }

  void Demon::post_integrate() {
    // Only run in simulation mode
    if (gflow->getRunMode()!=RunMode::SIM) return;

    // Compute door area
    RealType door_area = 0;
    if (size()>0) {
      // Update local variables
      update_local_ids(simData);
      // Look for bounding box.
      RealType mn = door_positions[0][1], mx = door_positions[0][1];
      for (const auto& v : door_positions) {        
        RealType y = v[1];
        if      (y<mn) mn = y;
        else if (mx<y) mx = y;
      }
      // Compute door area
      door_area = mx - mn;
    }

    // Check if the gas particles interact with one another.
    bool doesInteract = (forceMaster->getInteraction(0, 0)!=nullptr);
  
    RealType im = 0;
    int iter = 0, sz = simData->size();
    while (im==0 && iter<sz) {
      im = simData->Im(iter);
      ++iter;
    }

    // Record parameters
    parameters->addRecord("Vol", 0.5*gflow->getBounds().vol()); // This assumes the partitions have equal sizes.
    parameters->addRecord("Tau", tau);
    parameters->addRecord("Area", door_area);
    parameters->addRecord("Dimensions", sim_dimensions);
    parameters->addRecord("Interact", doesInteract ? 1. : 0.);
    if (im>0) parameters->addRecord("Mass", 1./im);
  }

  bool Demon::should_door_open(int nl, int nr, RealType el, RealType er) {
    return check_choice(nl, nr, el, er);
  }

  void Demon::setPartitionPosition(RealType x) {
    partition_position = x;
  }

  void Demon::setInteraction(DemonWall* w) {
    demon_interaction = w;
  }

  void Demon::setTau(RealType t) {
    if (t>0) tau = t;
  }

  inline void Demon::count_changes(int& nl, int& nr, RealType& el, RealType& er) {
    // Make sure we have a valid side entry array.
    if (side_entry<0) return;

    // Clear
    nl = nr = 0;
    el = er = 0.;
    Nl = Nr = 0;
    El = Er = 0.;

    // Get data from simdata.
    RealType **x = simData->X();
    RealType **v = simData->V();
    RealType *sg = simData->Sg();
    RealType *im = simData->Im();
    int *side = simData->IntegerData(side_entry);

    // Check all particles
    for (int i=0; i<simData->size(); ++i) {
      // Don't count walls.
      if (im[i]==0) continue;
      // Find current side
      int sd = get_side_literal(x[i]);
      // Calculate energy
      RealType ke = (1./im[i]) * 0.5 * sqr(v[i], sim_dimensions);
      // If the particle changed sides.
      if (sd!=side[i]) {
        // Left to right.
        if (sd==1) {
          ++nl;
          el += ke;
        }
        // Right to left.
        else if (sd==0) {
          ++nr;
          er += ke;
        }
        else;
      }

      // Count up particles and energy on each side.
      if (sd==0) {
        ++Nl;
        El += ke;
      }
      else {
        ++Nr;
        Er += ke;
      }
    }
  }

  inline void Demon::assign_side(bool nonliteral) {
    // Get data from simdata.
    RealType **x = simData->X();
    RealType *sg = simData->Sg();
    int *side = simData->IntegerData(side_entry);
    // Make sure side is a valid array.
    if (side==nullptr) return;
    // Check all particles
    for (int i=0; i<simData->size(); ++i) {
      // Find current side
      int sd;
      if (nonliteral) sd = get_side(x[i], side[i], sg[i]);
      else sd = get_side_literal(x[i]);
      // Write side data
      side[i] = sd;
    }
  }

  inline void Demon::checkpoint_data() {
    // Resize arrays if needbe
    if (checkpoint_x.size()!=simData->size()) {
      // Resize vectors
      checkpoint_x = vector<Vec>(simData->size(), Vec(sim_dimensions));
      checkpoint_v = vector<Vec>(simData->size(), Vec(sim_dimensions));
    }
    // Get arrays from simdata
    RealType **x = simData->X();
    RealType **v = simData->V();
    // Write simdata to arrays
    for (int i=0; i<simData->size(); ++i) {
      copyVec(x[i], checkpoint_x[i]);
      copyVec(v[i], checkpoint_v[i]);
    }

    // Checkpoint animation data.
    if (animate_object) {
      animate_last_recording = animate_object->getLastRecording();
      animate_last_size = animate_object->positions.size();
    }
  }

  inline void Demon::revert_to_last_checkpoint(bool construct) {
    // Check that sizes are the same
    if (simData->size()!=checkpoint_x.size()) throw false;

    // Get arrays from simdata
    RealType **x = simData->X();
    RealType **v = simData->V();

    // Write arrays to simdata
    for (int i=0; i<simData->size(); ++i) {
      copyVec(checkpoint_x[i], x[i]);
      copyVec(checkpoint_v[i], v[i]);
    }
    // Reset the time in gflow.
    gflow->setElapsedTime(last_check);

    // Force a rebuild of domains and forces
    if (construct) domain->construct();

    // Modify animation object.
    if (animate_object) {
      animate_object->setLastRecording(animate_last_recording);
      animate_object->positions.resize(animate_last_size);
      animate_object->timeStamps.resize(animate_last_size);
    }
  }

  inline void Demon::open_door() {
    // Where to put the particles.
    RealType pos = 2*gflow->getBounds().max[0];
    Vec p(sim_dimensions);
    p.set1(pos);
    // Update local ids.
    update_local_ids(simData);
    // Get data from simdata.
    RealType **x = simData->X();
    // Put particles far away
    for (int i=0; i<size(); ++i) {
      int id = at(i);
      // Put the particle back in place.
      copyVec(p, x[id]);
    }
    // Set door open
    door_open = true;
  }

  inline void Demon::close_door() {
    // Update local ids.
    update_local_ids(simData);
    // Get data from simdata.
    RealType **x = simData->X();
    // Put particles back
    for (int i=0; i<size(); ++i) {
      int id = at(i);
      // Put the particle back in place.
      copyVec(door_positions[i], x[id]);
    }
    // Set door open
    door_open = false;
    // Set demon wall to open
    if (demon_interaction) demon_interaction->turnOn();
    // Need to construct the domain. 
    // \todo Be able to add some particles to the domain without recreating the whole thing.
    domain->construct();
  }

  inline void Demon::add_data(int dn, RealType de) {
    dN.push_back(dn);
    dE.push_back(de);
  }

  inline int Demon::get_side(RealType *x, int current, RealType radius) {
    if (current==0) {
      if (partition_position + 2*radius <= x[0]) return 1;
      else return 0;
    }
    else {
      if (x[0] <= partition_position - 2*radius) return 0;
      else return 1;
    }
  }

  inline int Demon::get_side_literal(RealType *x) {
    return x[0]<=partition_position ? 0 : 1;
  }

  bool Demon::direction_demon(int nl, int nr, RealType el, RealType er) {
    return nr==0 && er==0.;
  }

  bool Demon::energy_demon(int nl, int nr, RealType el, RealType er) {
    return el-er>=0;
  }

  bool Demon::number_demon(int nl, int nr, RealType el, RealType er) {
    return nl-nr>=0;
  }

}