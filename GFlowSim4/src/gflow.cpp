#include "gflow.hpp"
#include "simdata.hpp"
#include "integrator.hpp"
#include "sectorization.hpp"
#include "neighbors.hpp"
#include "communicator.hpp"
#include "modifier.hpp"

namespace GFlowSimulation {

  GFlow::GFlow(int argc, char **argv) : running(false), fulfilled_time(0), iter(0) {
    simData = new SimData(this);
    sectorization = new Sectorization(this);
    neighbors = new Neighbors(this);
    communicator = new Communicator(this);
  }

  GFlow::~GFlow() {
    delete simData;
    delete sectorization;
    delete neighbors;
    delete communicator;
  }

  void GFlow::initialize() {
    simData->initialize();
    sectorization->initialize();
    neighbors->initialize();
    communicator->initialize();
  }

  void GFlow::run(RealType requested_time) {

    // Pre-integrate
    running = true;
    fulfilled_time = 0;
    iter = 0;

    // Do integration for the requested amount of time
    while (running) {
      // --- Pre-step
      integrator->pre_step();
      for (auto m : modifiers) m->pre_step();

      // --- Pre-exchange

      // --- Exchange particles (potentially) ---

      // --- Pre-neighbors

      // --- Do neighbor updates (potentially)

      // --- Pre-force
      integrator->pre_forces();
      for (auto m : modifiers) m->pre_forces();
      // Clip positions (if necessary)

      // --- Do forces

      // --- Post-forces
      integrator->post_forces();
      for (auto m : modifiers) m->post_forces();

      // Timer updates
      ++iter;

      // --- Post-step
      if (requested_time<=fulfilled_time) running = false;
      for (auto m : modifiers) m->post_step();
    }

    // Post-integrate

  }

}