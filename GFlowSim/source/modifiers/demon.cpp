#include "modifiers/demon.hpp"
// Other files
#include "base/interactionhandler.hpp"
#include "alldataobjects.hpp"
#include "base/datamaster.hpp"
#include "base/forcemaster.hpp"

#include "allinteractions.hpp"

using namespace GFlowSimulation;

Demon::Demon(GFlow *gflow)
    : Modifier(gflow), Group(gflow), startTime(0), tau(0.1) {
  // Set up data objects
  kineticL = make_shared<GraphObject>(gflow, "KineticL", "t", "KE");
  kineticR = make_shared<GraphObject>(gflow, "KineticR", "t", "KE");
  numberL = make_shared<GraphObject>(gflow, "NumberL", "t", "number");
  numberR = make_shared<GraphObject>(gflow, "NumberR", "t", "number");
  current_E = make_shared<GraphObject>(gflow, "CurrentE", "t", "dE");
  current_N = make_shared<GraphObject>(gflow, "CurrentN", "t", "dE");

  // Add data objects to gflow.
  gflow->addDataObject(kineticL);
  gflow->addDataObject(kineticR);
  gflow->addDataObject(numberL);
  gflow->addDataObject(numberR);
  gflow->addDataObject(current_E);
  gflow->addDataObject(current_N);

  // Create parameters object
  parameters = make_shared<Parameters>(gflow);
  gflow->addDataObject(parameters);

  // Default door choice function - direction demon.
  setDemon(0);
};

void Demon::pre_integrate() {
  // Only run in simulation mode
  if (gflow->getRunMode() != RunMode::SIM) {
    return;
  }

  // Get the side
  side_entry = simData->requestIntegerData("Side");
  // Entry must be positive. Otherwise, something is wrong.
  if (side_entry < 0) {
    throw false;
  }
  // Assign particles their side.
  assign_side();
  // Check point the initial data.
  checkpoint_data();
  // Set the last check time.
  last_check = gflow->getElapsedTime();

  // --- Record where door particles are

  // Set up vector of door particle positions.
  door_positions = vector<Vec>(size(), Vec(sim_dimensions));
  // Update local ids.
  update_local_ids();
  // Get data from simdata.
  auto x = simData->X();
  // Put particles back
  for (int i = 0; i < size(); ++i) {
    int id = at(i);
    // Put the particle back in place.
    copyVec(x(id), door_positions[i]);
  }

  // --- Look for an animation object, so we can modify it.
  for (auto obj : dataMaster->getDataObjects()) {
    auto ob = std::dynamic_pointer_cast<PositionData>(obj);
    if (ob) {
      animate_object = ob;
    }
  }
  // Reset
  opened_tau = closed_tau = 0;
  // Open the door
  open_door();
}

void Demon::pre_forces() {
  // Only run in simulation mode
  if (gflow->getRunMode() != RunMode::SIM) {
    return;
  }

  // Get the current time
  RealType time = gflow->getElapsedTime();

  // For counting changes
  int nl = 0, nr = 0;
  RealType el = 0, er = 0;

  // Check if enough time has gone by.
  if (time - last_check >= tau && startTime < time) {

    // If door was closed, then it was previously decided that the last time window should have been a closed window.
    // In this case, nl = nr = el = er = 0, and we can setup to check the next time window.
    // If the door was open, then we should count the changes and decide if we need to go back in time and keep the door
    // closed or not.
    if (door_open) {
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
        ++opened_tau;
      }
        // Otherwise, go back in time and keep the door closed.
      else {
        // Revert to last checkpoint. The particles' side data will be correct after this.
        revert_to_last_checkpoint(false); // Don't rebuild, since closing the door will rebuild.
        // Close the door.
        close_door();
        // Do not checkpoint data or assign side.
      }
    }
      // Door was closed.
    else {
      // Door was closed, so nl = nr = el = er = 0. Store this.
      add_data(0, 0);
      // Update El, Er, Nl, Nr.
      count_changes(nl, nr, el, er);
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
      ++closed_tau;
    }
  }
  else if (door_just_closed) {
    // Assign side the timestep after the door closes so the door can teleport particles to their proper sides.
    assign_side();
    // Set flag to false
    door_just_closed = false;
  }
}

void Demon::post_integrate() {
  // Only run in simulation mode
  if (gflow->getRunMode() != RunMode::SIM) {
    return;
  }

  // Compute door area
  RealType door_area = 0;
  if (size() > 0) {
    // Update local variables
    update_local_ids();
    // Look for bounding box.
    RealType mn = door_positions[0][1], mx = door_positions[0][1];
    for (const auto &v : door_positions) {
      RealType y = v[1];
      if (y < mn) {
        mn = y;
      }
      else if (mx < y) {
        mx = y;
      }
    }
    // Compute door area
    door_area = mx - mn;
  }

  // Check if the gas particles interact with one another.
  bool doesInteract = (forceMaster->getInteraction(0, 0) != nullptr);

  RealType im = 0;
  int iter = 0, sz = simData->size_owned();
  while (im == 0 && iter < sz) {
    im = simData->Im(iter);
    ++iter;
  }

  // Record parameters
  parameters->addRecord("Vol", 0.5 * gflow->getBounds().vol()); // This assumes the partitions have equal sizes.
  parameters->addRecord("Tau", tau);
  parameters->addRecord("Area", door_area);
  parameters->addRecord("Dimensions", sim_dimensions);
  parameters->addRecord("Interact", doesInteract ? 1. : 0.);
  parameters->addRecord("Type", demon_type);
  if (im > 0) {
    parameters->addRecord("Mass", 1. / im);
  }
}

bool Demon::should_door_open(int nl, int nr, RealType el, RealType er) {
  return check_choice(nl, nr, el, er);
}

void Demon::setPartitionPosition(RealType x) {
  partition_position = x;
}

void Demon::setInteraction(shared_ptr<DemonWall> w) {
  demon_interaction = w;
}

void Demon::setStartTime(RealType t) {
  startTime = t;
}

void Demon::setTau(RealType t) {
  if (t > 0) {
    tau = t;
  }
}

void Demon::setDemon(int choice) {
  if (0 <= choice && choice <= 2) {
    demon_type = choice;
    if (choice == 0) {
      check_choice = &direction_demon;
    }
    else if (choice == 1) {
      check_choice = &energy_demon;
    }
    else if (choice == 2) {
      check_choice = &number_demon;
    }
  }
}

inline void Demon::count_changes(int &nl, int &nr, RealType &el, RealType &er) {
  // Make sure we have a valid side entry array.
  if (side_entry < 0) {
    return;
  }

  // Clear
  nl = nr = 0;
  el = er = 0.;
  Nl = Nr = 0;
  El = Er = 0.;

  // Get data from simdata.
  auto x = simData->X();
  auto v = simData->V();
  auto im = simData->Im();
  auto side = simData->IntegerData(side_entry);

  // Check all particles
  for (int i = 0; i < simData->size_owned(); ++i) {
    // Don't count walls.
    if (im(i) == 0) {
      continue;
    }
    // Find current side
    int sd = get_side(x(i));
    // Calculate energy
    RealType ke = (1. / im(i)) * 0.5 * sqr(v(i), sim_dimensions);
    // If the particle changed sides.
    if (sd != side(i)) {
      // Left to right.
      if (sd == 1) {
        ++nl;
        el += ke;
      }
        // Right to left.
      else if (sd == 0) {
        ++nr;
        er += ke;
      }
    }
    // Count up particles and energy on each side.
    if (sd == 0) {
      ++Nl;
      El += ke;
    }
    else {
      ++Nr;
      Er += ke;
    }
  }
}

inline void Demon::assign_side() {
  // Get data from simdata.
  auto x = simData->X();
  auto side = simData->IntegerData(side_entry);
  // Make sure side is a valid array.
  if (side.isnull()) {
    return;
  }
  // Check all particles
  for (int i = 0; i < simData->size_owned(); ++i) {
    // Find current side
    int sd = get_side(x(i));
    // Write side data
    side(i) = sd;
  }
}

inline void Demon::checkpoint_data() {
  // Resize arrays if needbe
  if (checkpoint_x.size() != simData->size_owned()) {
    // Resize vectors
    checkpoint_x = vector<Vec>(simData->size_owned(), Vec(sim_dimensions));
    checkpoint_v = vector<Vec>(simData->size_owned(), Vec(sim_dimensions));
  }
  // Get arrays from simdata
  auto x = simData->X();
  auto v = simData->V();
  // Write simdata to arrays
  for (int i = 0; i < simData->size_owned(); ++i) {
    copyVec(x(i), checkpoint_x[i]);
    copyVec(v(i), checkpoint_v[i]);
  }

  // Checkpoint animation data.
  if (animate_object) {
    animate_last_recording = animate_object->getLastRecording();
    animate_last_size = animate_object->positions.size();
  }
}

inline void Demon::revert_to_last_checkpoint(bool construct) {
  // Check that sizes are the same
  if (simData->size_owned() != checkpoint_x.size()) {
    throw false;
  }

  // Get arrays from simdata
  auto x = simData->X();
  auto v = simData->V();

  // Write arrays to simdata
  for (int i = 0; i < simData->size_owned(); ++i) {
    copyVec(checkpoint_x[i], x(i));
    copyVec(checkpoint_v[i], v(i));
  }
  // Reset the time in gflow.
  gflow->setElapsedTime(last_check);

  // Force a rebuild of handler and forces
  if (construct) {
    handler->construct();
  }

  // Modify animation object.
  if (animate_object) {
    animate_object->setLastRecording(animate_last_recording);
    animate_object->positions.resize(animate_last_size);
    animate_object->timeStamps.resize(animate_last_size);
  }
}

inline void Demon::open_door() {
  // Where to put the particles.
  RealType pos = 2 * gflow->getBounds().max[0];
  Vec p(sim_dimensions);
  p.set1(pos);
  // Update local ids.
  update_local_ids();
  // Get data from simdata.
  auto x = simData->X();
  // Put particles far away
  for (int i = 0; i < size(); ++i) {
    int id = at(i);
    // Put the particle back in place.
    copyVec(p, x(id));
  }
  // Set door open
  door_open = true;
}

inline void Demon::close_door() {
  // Update local ids.
  update_local_ids();
  // Get data from simdata.
  auto x = simData->X();
  // Put particles back
  for (int i = 0; i < size(); ++i) {
    int id = at(i);
    // Put the particle back in place.
    copyVec(door_positions[i], x(id));
  }
  // Set door open
  door_open = false;
  // Set demon wall to open
  if (demon_interaction) {
    demon_interaction->turnOn();
  }
  // Need to construct the handler.
  // \todo Be able to add some particles to the handler without recreating the whole thing.
  handler->construct();
  // Set the flag
  door_just_closed = true;
}

inline void Demon::add_data(int dn, RealType de) {
  // Get the current time
  RealType time = gflow->getElapsedTime();
  current_E->addEntry(time, de);
  current_N->addEntry(time, dn);
}

inline int Demon::get_side(RealType *x) {
  return x[0] <= partition_position ? 0 : 1;
}

bool Demon::direction_demon(int nl, int nr, RealType el, RealType er) {
  return nr == 0 && er == 0.;
}

bool Demon::energy_demon(int nl, int nr, RealType el, RealType er) {
  return el >= er;
}

bool Demon::number_demon(int nl, int nr, RealType el, RealType er) {
  return nl >= nr;
}
