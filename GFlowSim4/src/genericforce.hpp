#ifndef __GENERIC_FORCE_HPP__GFLOW__
#define __GENERIC_FORCE_HPP__GFLOW__

#include "interaction.hpp"

namespace GFlowSimulation {

  //! @todo Implement this
  template<int NTerms> class GenericForce : public Interaction {
    // Constructor
    GenericForce(GFlow *);
  };

}
#endif // __GENERIC_FORCE_HPP__GFLOW__