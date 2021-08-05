#include "interactions/interaction-choice.hpp"
// Other files.
#include "interactions/demon_wall.hpp"

using namespace GFlowSimulation;

shared_ptr<Interaction> InteractionChoice::choose(GFlow *gflow,
                                                  const string &token,
                                                  int sim_dimensions) {
  string err = " Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in";
  if (token == NoneToken) {
    return nullptr;
  }
  else if (token == HardSphereToken) {
    if (sim_dimensions == 1) {
      return make_shared<HardSphereVLP<1> >(gflow);
    }
    else if (sim_dimensions == 2) {
      return make_shared<HardSphereVLP<2> >(gflow);
    }
    else if (sim_dimensions == 3) {
      return make_shared<HardSphereVLP<3> >(gflow);
    }
    else if (sim_dimensions == 4) {
      return make_shared<HardSphereVLP<4> >(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
  }
  else if (token == HardSphereDsToken) {
    if (sim_dimensions == 1) {
      return make_shared<HardSphereDsVLP<1> >(gflow);
    }
    else if (sim_dimensions == 2) {
      return make_shared<HardSphereDsVLP<2> >(gflow);
    }
    else if (sim_dimensions == 3) {
      return make_shared<HardSphereDsVLP<3> >(gflow);
    }
    else if (sim_dimensions == 4) {
      return make_shared<HardSphereDsVLP<4> >(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
  }
  else if (token == LennardJonesToken) {
    if (sim_dimensions == 1) {
      return make_shared<LennardJonesVLP<1> >(gflow);
    }
    else if (sim_dimensions == 2) {
      return make_shared<LennardJonesVLP<2> >(gflow);
    }
    else if (sim_dimensions == 3) {
      return make_shared<LennardJonesVLP<3> >(gflow);
    }
    else if (sim_dimensions == 4) {
      return make_shared<LennardJonesVLP<4> >(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
  }
  else if (token == CoulombToken) {
    if (sim_dimensions == 1) {
      return make_shared<CoulombVLP<1> >(gflow);
    }
    else if (sim_dimensions == 2) {
      return make_shared<CoulombVLP<2> >(gflow);
    }
    else if (sim_dimensions == 3) {
      return make_shared<CoulombVLP<3> >(gflow);
    }
    else if (sim_dimensions == 4) {
      return make_shared<CoulombVLP<4> >(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
  }
  else if (token == BuckinghamToken) {
    if (sim_dimensions == 1) {
      return make_shared<BuckinghamVLP<1> >(gflow);
    }
    else if (sim_dimensions == 2) {
      return make_shared<BuckinghamVLP<2> >(gflow);
    }
    else if (sim_dimensions == 3) {
      return make_shared<BuckinghamVLP<3> >(gflow);
    }
    else if (sim_dimensions == 4) {
      return make_shared<BuckinghamVLP<4> >(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
  }
  else if (token == HertzToken) {
    if (sim_dimensions == 1) {
      return make_shared<HertzVLP<1> >(gflow);
    }
    else if (sim_dimensions == 2) {
      return make_shared<HertzVLP<2> >(gflow);
    }
    else if (sim_dimensions == 3) {
      return make_shared<HertzVLP<3> >(gflow);
    }
    else if (sim_dimensions == 4) {
      return make_shared<HertzVLP<4> >(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
  }
  else if (token == HookeToken) {
    if (sim_dimensions == 1) {
      return make_shared<HookeVLP<1> >(gflow);
    }
    else if (sim_dimensions == 2) {
      return make_shared<HookeVLP<2> >(gflow);
    }
    else if (sim_dimensions == 3) {
      return make_shared<HookeVLP<3> >(gflow);
    }
    else if (sim_dimensions == 4) {
      return make_shared<HookeVLP<4> >(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
  }

    // Non-generic types.
  else if (token == DemonWallToken) {
    if (sim_dimensions == 2) {
      return make_shared<DemonWall>(gflow);
    }
    else {
      throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
  }
  else {
    throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
  }
  // This point is never reached
}
