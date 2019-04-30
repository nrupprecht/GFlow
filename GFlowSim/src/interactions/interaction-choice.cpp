#include "interaction-choice.hpp"

#include "demon_wall.hpp"

namespace GFlowSimulation {

  Interaction* InteractionChoice::choose(GFlow *gflow, const string &token, int sim_dimensions) {
    if (token==NoneToken) return nullptr;
    else if (token==HardSphereToken) {
      if (sim_dimensions==2)
        return new HardSphere_VerletPairs_2d(gflow);
      else if (sim_dimensions==3)
        return new HardSphere_VerletPairs_3d(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else if (token==HardSphereDsToken) {
      if (sim_dimensions==2)
        return new HardSphereDs_VerletPairs_2d(gflow);
      else if (sim_dimensions==3)
        return new HardSphereDs_VerletPairs_3d(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else if (token==LennardJonesToken) {
      if (sim_dimensions==2)
        return new LennardJones_VerletPairs_2d(gflow);
      else if (sim_dimensions==3)
        return new LennardJones_VerletPairs_3d(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else if (token==BuckinghamToken) {
      if (sim_dimensions==2)
        return new Buckingham_VerletPairs_2d(gflow);
      else if (sim_dimensions==3)
        return new Buckingham_VerletPairs_3d(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else if (token==DetectorToken) {
      if (sim_dimensions==2)
        return new Detector_2d(gflow);
      if (sim_dimensions==3)
        return new Detector_3d(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else if (token==HardSphereReflectingToken) {
      if (sim_dimensions==2)
        return new HardSphere_Reflecting_2d(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else if (token==DemonWallToken) {
      if (sim_dimensions==2)
        return new DemonWall(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    // This point is never reached
  }

}