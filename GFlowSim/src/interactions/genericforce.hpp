#ifndef __GENERIC_FORCE_HPP__GFLOW__
#define __GENERIC_FORCE_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  //! @todo Implement this
  template<int NTerms> class GenericForce : public Interaction {
    // Constructor
    GenericForce(GFlow *);

    //! @brief Initialize the force, check if all the special data (dataF, dataI) the force needs exists, make
    //! sure parameter packs are up to date.
    virtual void initialize() override;
  };

}
#endif // __GENERIC_FORCE_HPP__GFLOW__