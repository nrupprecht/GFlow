#ifndef __STRIPE_X_HPP__GFLOW__
#define __STRIPE_X_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class StripeX : public Modifier {
  public:
    //! \brief Default constructor.
    StripeX(GFlow*);

    //! \brief Have simdata set up a special data array for each particle.
    virtual void pre_integrate();

    //! \brief Set the special data array as necessary.
    virtual void post_forces();

  private:

    //! \brief Which entry simdata saves the stripe data in.
    int entry = -1;

    //! \brief How wide the window is that we consider the setting part.
    RealType window;

    // So we don't apply the force all the time - this is taxing because of all the random numbers
    // we need to generate
    RealType lastUpdate, updateDelay;
  };

}
#endif // __STRIPE_X_HPP__GFLOW__