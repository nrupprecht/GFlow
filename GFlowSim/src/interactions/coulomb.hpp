#ifndef __COULOMB_HPP__GFLOW__
#define __COULOMB_HPP__GFLOW__

#include "../base/interaction.hpp"
#include "interaction__2d.hpp"

namespace GFlowSimulation {

  class Coulomb2d : public Interaction2d {
  public:
    //! \brief Default constructor.
    Coulomb2d(GFlow*);

    //! \brief Set the cutoff, and change the potential energy shift.
    void setCutoff(RealType);

    //! \brief Suggests a safe timescale given the minimum mass of a particle that has this interaction.
    RealType suggest_timescale(RealType mass) const override {
      return -1; // \todo Calculate this.
    }

  private:
    //! \brief The interaction kernel.
    virtual void kernel(int, int, RealType, RealType, RealType, RealType*, RealType**) const override;

    //! \brief The strength of the coulomb repulsion.
    RealType repulsion;

    //! \brief How much we have to shift the potential energy by to take into account the fact that there
    //! is a force cutoff.
    RealType potential_energy_shift;
  };

}
#endif // __COULOMB_HPP__GFLOW__