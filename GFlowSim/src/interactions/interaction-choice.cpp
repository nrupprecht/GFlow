#include "interaction-choice.hpp"

#include "demon_wall.hpp"
    
namespace GFlowSimulation {
  
  Interaction* InteractionChoice::choose(GFlow *gflow, const string &token, int sim_dimensions) {
    string err = " Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in";
    if (token==NoneToken) return nullptr;
    else if (token==HardSphereToken) {
      if      (sim_dimensions==1) return new HardSphereVLP<1>(gflow);
      else if (sim_dimensions==2) return new HardSphereVLP<2>(gflow);
      else if (sim_dimensions==3) return new HardSphereVLP<3>(gflow);
      else if (sim_dimensions==4) return new HardSphereVLP<4>(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
    else if (token==HardSphereDsToken) {
      if      (sim_dimensions==1) return new HardSphereDsVLP<1>(gflow);
      else if (sim_dimensions==2) return new HardSphereDsVLP<2>(gflow);
      else if (sim_dimensions==3) return new HardSphereDsVLP<3>(gflow);
      else if (sim_dimensions==4) return new HardSphereDsVLP<4>(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
    else if (token==LennardJonesToken) {
      if      (sim_dimensions==1) return new LennardJonesVLP<1>(gflow);
      else if (sim_dimensions==2) return new LennardJonesVLP<2>(gflow);
      else if (sim_dimensions==3) return new LennardJonesVLP<3>(gflow);
      else if (sim_dimensions==4) return new LennardJonesVLP<4>(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
    else if (token==CoulombToken) {
      if      (sim_dimensions==1) return new CoulombVLP<1>(gflow);
      else if (sim_dimensions==2) return new CoulombVLP<2>(gflow); 
      else if (sim_dimensions==3) return new CoulombVLP<3>(gflow);
      else if (sim_dimensions==4) return new CoulombVLP<4>(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
    else if (token==BuckinghamToken) {
      if      (sim_dimensions==1) return new BuckinghamVLP<1>(gflow);
      else if (sim_dimensions==2) return new BuckinghamVLP<2>(gflow);
      else if (sim_dimensions==3) return new BuckinghamVLP<3>(gflow);
      else if (sim_dimensions==4) return new BuckinghamVLP<4>(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
    else if (token==HertzToken) {
      if      (sim_dimensions==1) return new HertzVLP<1>(gflow);
      else if (sim_dimensions==2) return new HertzVLP<2>(gflow);
      else if (sim_dimensions==3) return new HertzVLP<3>(gflow);
      else if (sim_dimensions==4) return new HertzVLP<4>(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }
    else if (token==HookeToken) {
      if      (sim_dimensions==1) return new HookeVLP<1>(gflow);  
      else if (sim_dimensions==2) return new HookeVLP<2>(gflow);
      else if (sim_dimensions==3) return new HookeVLP<3>(gflow);
      else if (sim_dimensions==4) return new HookeVLP<4>(gflow);
      // Nothing per-say invalid about requesting this force in this dimension, it just isn't coded in.
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions) + "." + err);
    }

    // Non-generic types.
    else if (token==DemonWallToken) {
      if (sim_dimensions==2)
        return new DemonWall(gflow);
      else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    }
    else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
    // This point is never reached
  }

}
