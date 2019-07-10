#ifndef __BUCKINGHAM__2D_HPP__GFLOW__
#define __BUCKINGHAM__2D_HPP__GFLOW__

#include "buckingham.hpp"

namespace GFlowSimulation {

  class Buckingham_2d : public Buckingham {
  public:
    //! \brief Default constructor.
    Buckingham_2d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}
#endif // __BUCKINGHAM__2D_HPP__GFLOW__