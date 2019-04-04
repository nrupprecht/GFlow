#ifndef __INTERACTION_CHOICE_HPP__GFLOW__
#define __INTERACTION_CHOICE_HPP__GFLOW__

#include "../allinteractions.hpp"

namespace GFlowSimulation {
  namespace InteractionChoice {

    const string NoneToken = "None";
    const string HardSphereToken = "HardSphere";
    const string HardSphereDsToken = "HardSphereDs";
    const string LennardJonesToken = "LennardJones";
    const string BuckinghamToken = "Buckingham";
    const string DetectorToken = "Detector";

    /**
    *  \brief Error class for trying to pick an interaction that doesn't exist.
    *
    */
    class InvalidInteraction : public Exception {
    public:
      //! \brief Default constructor
      InvalidInteraction() : Exception() {};
      //! \brief Message constructor
      InvalidInteraction(const string& message) : Exception(message) {};
    };

    //! \brief Choose an interaction given a token and a dimensionality.
    inline Interaction* choose(GFlow *gflow, const string &token, int sim_dimensions) {
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
        if (sim_dimensions==3)
          return new Detector(gflow);
        else throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
      }
      else {
        throw InvalidInteraction(token + ", " + toStr(sim_dimensions));
      }
      // This point is never reached
    }

  }
}
#endif // __INTERACTION_CHOICE_HPP__GFLOW__