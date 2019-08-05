#ifndef __HERTZ_FORCE_HPP__GFLOW__
#define __HERTZ_FORCE_HPP__GFLOW__

namespace GFlowSimulation {

  //! \brief Hertzian particle interaction.
  //!
  //! Normal force:     \vec{F} = sqrt(r/D)*(K_n * r * \hat{n} - m/2 gamma_n * \vec{v}_n)
  //! Tangential force:
  class HertzForce : public Interaction {
  public:

  protected:
    //! \brief Normal elastic constant.
    RealType K_n;
    //! \brief Tangential elastic constant.
    RealType K_t;

    //! \brief Normal viscoelastic constant.
    RealType gamma_n;
    //! \brief Tangential viscoelastic constant.
    RealType gamma_t;
  };

}
#endif // __HERTZ_FORCE_HPP__GFLOW__