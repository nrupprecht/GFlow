#ifndef __INTERACTION_CHOICE_HPP__GFLOW__
#define __INTERACTION_CHOICE_HPP__GFLOW__

#include "../allinteractions.hpp"

namespace GFlowSimulation {

  const string NoneToken = "None";
  const string HardSphereToken = "HardSphere";
  const string HardSphereDsToken = "HardSphereDs";
  const string HardSphereCfToken = "HardSphereCf";
  const string LennardJonesToken = "LennardJones";
  const string BuckinghamToken = "Buckingham";
  const string HardSphereReflectingToken = "HardSphereReflecting";
  const string DemonWallToken = "DemonWall";
  const string CoulombToken = "Coulomb";
  const string HertzToken = "Hertz";
  const string HookeToken = "Hooke";

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

  class InteractionChoice {
  public:
    //! \brief Choose an interaction given a token and a dimensionality.
    static Interaction* choose(GFlow *gflow, const string &token, int sim_dimensions);
  };

}
#endif // __INTERACTION_CHOICE_HPP__GFLOW__