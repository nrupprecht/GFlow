#ifndef __HARD_SPHERE_DISSIPATIVE_HPP__GFLOW__
#define __HARD_SPHERE_DISSIPATIVE_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  /** 
  *  @brief A very general hard sphere force where all particles have their own repulsion, coefficient of friction, and dissipation.
  *
  *  The force between particles is proportional to their overlap, (r1 + r2 - distance), with constant of proportionality
  *  0.5(repulsion1 + repulsion2). There is a frictional force between the spheres (a tangential force), and a dissipative force 
  *  that slows down spheres moving towards one another.
  *
  */
  class GeneralHardSphere : public Interaction {
  public:
    //! @brief Constructor
    GeneralHardSphere(GFlow *);

    //! @brief Initialize the force, check if all the special data (dataF, dataI) the force needs exists, make
    //! sure parameter packs are up to date.
    virtual void initialize() override;

    //! @param[in] normal
    //! @param[in] distance
    //! @param[in] id1
    //! @param[in] id2
    //! @param[in] simData
    //! @param[in] param_pack A parameter pack, passed in from force. Contains characteristic 
    //! constants of the force, and extra data the force needs.
    //! @param[in,out] data_pack Data to be updated by the function.
    static void force(RealType*, const RealType, const int, const int, SimData*, const RealType*, RealType*);

  private:
    //! @brief Get the places from simdata using the labels.
    void getPlaces();

    string repulsionLabel, dissipationLabel, frictionLabel;
    int repulsionPlace, dissipationPlace, frictionPlace;

  };

}
#endif // __HARD_SPHERE_DISSIPATIVE_HPP__GFLOW__